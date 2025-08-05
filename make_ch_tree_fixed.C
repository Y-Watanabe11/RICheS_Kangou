// 使い方（対話入力）:  root -l -q 'make_ch_tree_fixed.C()'
// 直接指定:            root -l -q 'make_ch_tree_fixed.C("path/to/data.csv","out.root",true)'

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TObjString.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <cctype>

static std::string trim(const std::string& s){
    auto a = s.begin(), b = s.end();
    while (a!=b && isspace((unsigned char)*a)) ++a;
    while (b!=a && isspace((unsigned char)*(b-1))) --b;
    std::string t(a,b);
    if (!t.empty() && (t.front()=='"' || t.front()=='\'')) t.erase(t.begin());
    if (!t.empty() && (t.back()=='"'  || t.back()=='\'' )) t.pop_back();
    return t;
}
static inline std::string trim_field(const std::string& s){
    size_t a = 0, b = s.size();
    while (a<b && isspace((unsigned char)s[a])) ++a;
    while (b>a && isspace((unsigned char)s[b-1])) --b;
    return s.substr(a, b-a);
}
static std::string default_rootname_from_csv(const std::string& csv){
    auto p = csv.find_last_of("/\\");
    std::string base = (p==std::string::npos) ? csv : csv.substr(p+1);
    auto q = base.find_last_of('.');
    if (q!=std::string::npos) base = base.substr(0,q);
    return base + ".root";
}

// ホワイトスペース分割を強制（行内の , ; \t は空白に置換）
static void tokenize_ws_numbers(std::string line, std::vector<double>& out_vals){
    // UTF-8 BOM除去
    if (line.size()>=3 &&
        (unsigned char)line[0]==0xEF && (unsigned char)line[1]==0xBB && (unsigned char)line[2]==0xBF){
        line.erase(0,3);
    }
    for (char& c : line){
        if (c==',' || c==';' || c=='\t') c = ' ';
    }
    std::istringstream iss(line);
    std::string tok;
    out_vals.clear();
    while (iss >> tok){
        tok = trim_field(tok);
        try { out_vals.push_back(std::stod(tok)); }
        catch (...) { out_vals.push_back(0.0); }
    }
}

void make_ch_tree_fixed(const char* csvfile_cstr = "",
                        const char* outroot_cstr = "",
                        bool flipY = true)  // 見た目上 y を上向きに
{
    // 列名を固定（ヘッダー無し前提）
    std::vector<std::string> COLS;
    COLS.reserve(5+64);
    COLS.emplace_back("Dummy_flag");
    COLS.emplace_back("OUT_flag");
    COLS.emplace_back("Event_Flag");
    COLS.emplace_back("count");
    COLS.emplace_back("time");
    for (int i=1;i<=64;++i) COLS.emplace_back(Form("ch%d", i));
    const int N = (int)COLS.size(); // 69

    // 入力取得
    std::string csv = trim(csvfile_cstr ? std::string(csvfile_cstr) : std::string());
    if (csv.empty()) {
        std::cout << "CSVファイルのパスを入力してください（ドラッグ&ドロップ可）: ";
        std::getline(std::cin, csv);
        csv = trim(csv);
    }
    if (csv.empty()) { std::cerr << "[err] CSVパスが空です\n"; return; }

    std::string outroot = trim(outroot_cstr ? std::string(outroot_cstr) : std::string());
    if (outroot.empty()) {
        std::string def = default_rootname_from_csv(csv);
        std::cout << "出力ROOTファイル名（既定: " << def << "）: ";
        std::string in; std::getline(std::cin, in);
        in = trim(in);
        outroot = in.empty() ? def : in;
    }

    // 入出力
    std::ifstream fin(csv);
    if (!fin) { std::cerr << "[err] open fail: " << csv << "\n"; return; }

    std::unique_ptr<TFile> fout(TFile::Open(outroot.c_str(), "RECREATE"));
    if (!fout || fout->IsZombie()) { std::cerr << "[err] cannot create: " << outroot << "\n"; return; }

    // メタ：列名メモ
    {
        TDirectory* meta = fout->mkdir("meta"); meta->cd();
        TString note = "columns=Dummy_flag,OUT_flag,Event_Flag,count,time";
        for (int i=1;i<=64;++i) note += Form(",ch%d", i);
        (new TObjString(note))->Write("columns_note");
        fout->cd();
    }

    // TTree
    TTree* tree = new TTree("ch_eventtree", "channel-wise values per csv row");
    Long64_t event = 0;                      // 0行目→event=0
    tree->Branch("event", &event);

    std::vector<double> vbuf(N, 0.0);
    for (int j=0;j<N;++j){
        tree->Branch(COLS[j].c_str(), &vbuf[j], Form("%s/D", COLS[j].c_str()));
    }

    // ヒスト保存用ディレクトリ
    TDirectory* dh = fout->mkdir("hist");

    // CSV読み込み（0行目からイベントとして使用）
    const int nx=8, ny=8;
    const int start0 = 5;        // ch1 の開始列(0始まり) = 6列目
    const int need   = nx*ny;    // 64

    std::string line;
    std::vector<double> vals;
    Long64_t nSavedHist=0;
    Long64_t lineNo=0;

    while (std::getline(fin, line)) {
        // 空行スキップ（event番号は進めない）
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

        tokenize_ws_numbers(line, vals);
        ++lineNo;

        // デバッグ（最初の数行だけ）
        if (lineNo <= 3) {
            std::cout << "[info] line " << lineNo << " tokens=" << vals.size() << std::endl;
        }

        // vbufに格納（不足は0、超過は切捨て）
        for (int j=0;j<N;++j) vbuf[j] = (j<(int)vals.size()) ? vals[j] : 0.0;

        tree->Fill();                  // 0行目が event=0

        // 8×8 ヒスト保存（ch1..ch64）
        if ((int)vals.size() >= start0 + need) {
            dh->cd();
            TH2D* h = new TH2D(Form("histtree_%06lld", event),
                               Form("Event %lld;X bin;Y bin", event),
                               nx, 0.5, nx+0.5, ny, 0.5, ny+0.5);
            for (int iy=0; iy<ny; ++iy) {
                for (int ix=0; ix<nx; ++ix) {
                    double v = vbuf[start0 + iy*nx + ix];
                    int ybin = flipY ? (ny - iy) : (iy + 1);
                    h->SetBinContent(ix+1, ybin, v);
                }
            }
            h->Write(); fout->cd(); ++nSavedHist;
        } else {
            std::cerr << "[warn] line " << lineNo
                      << " has only " << vals.size()
                      << " tokens (< " << (start0+need) << "), hist skipped.\n";
        }
        ++event;                        // 次の行は event+1
    }

    fout->cd(); tree->Write(); fout->Write();
    std::cout << "[done] wrote " << outroot
              << "  entries=" << tree->GetEntries()
              << "  hist_saved=" << nSavedHist << std::endl;
}
