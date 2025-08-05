// 実行例（対話入力モード）:
//   root -l -q 'draw_event_hist_fixed.C()'
// 直接指定（自動保存ON・名前既定・PNG保存なし）:
//   root -l -q 'draw_event_hist_fixed.C("out.root",123,true,-1,-1,"",true,"hist","", "g", 1.2)'

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <iostream>
#include <memory>
#include <string>

static std::string trim2(const std::string& s){
    auto a = s.begin(), b = s.end();
    while (a!=b && isspace((unsigned char)*a)) ++a;
    while (b!=a && isspace((unsigned char)*(b-1))) --b;
    std::string t(a,b);
    if (!t.empty() && (t.front()=='"' || t.front()=='\'')) t.erase(t.begin());
    if (!t.empty() && (t.back()=='"'  || t.back()=='\'' )) t.pop_back();
    return t;
}

static TH2D* build_from_tree(TTree* tree, Long64_t ev, bool flipY) {
    const int nx=8, ny=8;
    if (!tree || ev<0 || ev>=tree->GetEntries()) return nullptr;
    tree->GetEntry(ev);
    auto* h = new TH2D(Form("histtree_%06lld_tmp", ev),
                       Form("Event %lld;X bin;Y bin", ev),
                       nx, 0.5, nx+0.5, ny, 0.5, ny+0.5);
    for (int iy=0; iy<ny; ++iy){
        for (int ix=0; ix<nx; ++ix){
            int ch = iy*nx + ix + 1; // ch1..ch64
            TLeaf* L = tree->GetLeaf(Form("ch%d", ch));
            double v = L ? L->GetValue() : 0.0;
            int ybin = flipY ? (ny - iy) : (iy + 1);
            h->SetBinContent(ix+1, ybin, v);
        }
    }
    return h;
}

static TDirectory* ensure_dir(TFile* f, const char* dname){
    if (!f) return nullptr;
    TDirectory* d = (TDirectory*)f->Get(dname);
    if (!d) d = f->mkdir(dname);
    return d;
}

void draw_event_hist_fixed(const char* rootfile_cstr = "",
                           Long64_t event = -1,
                           bool flipY = true,
                           double zmin = -1, double zmax = -1,
                           const char* save_png_cstr = "",
                           bool save_to_root = false,           // 自動保存するか
                           const char* hist_dir_cstr = "hist",  // 保存先ディレクトリ名
                           const char* hist_name_cstr = "",     // 保存名（空なら histtree_%06lld）
                           const char* paintTextFormat = "g",   // ← %なし ("g", ".0f", "4.1f" 等)
                           double textMarkerSize = 1.2)         // 数値テキストの大きさ
{
    // 対話入力（空・未指定のとき）
    std::string rootfile = trim2(rootfile_cstr ? std::string(rootfile_cstr) : std::string());
    if (rootfile.empty()){
        std::cout << "ROOTファイルのパスを入力してください: ";
        std::getline(std::cin, rootfile);
        rootfile = trim2(rootfile);
        if (rootfile.empty()){ std::cerr << "[err] ROOTファイル未指定\n"; return; }
    }
    if (event < 0){
        std::cout << "表示する event 番号を入力してください（0始まり）: ";
        std::string s; std::getline(std::cin, s);
        try { event = std::stoll(s); } catch (...) { std::cerr << "[err] 数値で入力してください\n"; return; }
    }

    // 読み書き両用で開く（UPDATE）
    std::unique_ptr<TFile> f(TFile::Open(rootfile.c_str(), "UPDATE"));
    if (!f || f->IsZombie()){ std::cerr << "[err] open fail: " << rootfile << "\n"; return; }

    // 既存ヒスト（hist/以下）を優先取得
    TH2D* h = nullptr;
    if (auto* dh = (TDirectory*)f->Get("hist")){
        h = (TH2D*)dh->Get(Form("histtree_%06lld", event));
    }

    // 無ければ ch_eventtree から生成
    std::unique_ptr<TH2D> holder;
    if (!h){
        TTree* tree = (TTree*)f->Get("ch_eventtree");
        if (!tree){ std::cerr << "[err] ch_eventtree not found\n"; return; }
        holder.reset(build_from_tree(tree, event, flipY));
        if (!holder){ std::cerr << "[err] cannot build hist for event " << event << "\n"; return; }
        h = holder.get();
    }

    // 描画設定（色 + 数値）
    gStyle->SetNumberContours(50);
    gStyle->SetPaintTextFormat(paintTextFormat); // 例: "g", ".0f", "4.1f"
    gStyle->SetOptStat(0);                       // 統計箱オフ（必要なら外す）
    h->SetMarkerSize(textMarkerSize);

    TCanvas* c = (TCanvas*)gROOT->FindObject("c_event");
    if (!c) c = new TCanvas("c_event","event view",1000,850);
    c->cd();

    if (zmin>=0) h->SetMinimum(zmin);
    if (zmax>=0) h->SetMaximum(zmax);

    // カラーマップ + 数値を同時描画
    h->Draw("COLZ TEXT");

    // タイトル（Event番号）
    TLatex lat; lat.SetNDC(true); lat.SetTextSize(0.045);
    lat.DrawLatex(0.14, 0.945, Form("#bf{Event %lld}", event));
    c->Update();

    // PNG保存（対話）
    std::string save_png = trim2(save_png_cstr ? std::string(save_png_cstr) : std::string());
    if (save_png.empty()){
        std::cout << "PNG保存ファイル名（空で保存しない）: ";
        std::getline(std::cin, save_png);
        save_png = trim2(save_png);
    }
    if (!save_png.empty()){
        c->SaveAs(save_png.c_str());
        std::cout << "[saved] " << save_png << "\n";
    }

    // ROOTへ保存（インタラクティブ or 自動）
    std::string hist_dir  = trim2(hist_dir_cstr ? std::string(hist_dir_cstr) : std::string("hist"));
    std::string hist_name = trim2(hist_name_cstr ? std::string(hist_name_cstr) : std::string());
    if (!save_to_root){
        std::cout << "このヒストをROOTファイルに保存しますか？ [y/N]: ";
        std::string yn; std::getline(std::cin, yn);
        if (!yn.empty() && (yn=="y" || yn=="Y")) save_to_root = true;
    }
    if (save_to_root){
        if (hist_name.empty()) hist_name = Form("histtree_%06lld", event);
        TDirectory* dh = ensure_dir(f.get(), hist_dir.c_str());
        if (!dh){ std::cerr << "[err] cannot create/find dir: " << hist_dir << "\n"; return; }
        dh->cd();

        // 既存があれば上書き
        TH2D* hexist = (TH2D*)dh->Get(hist_name.c_str());
        if (hexist) dh->Delete((hist_name+";*").c_str());

        TH2D* hsave = (TH2D*)h->Clone(hist_name.c_str());
        hsave->SetTitle(Form("Event %lld;X bin;Y bin", event));
        hsave->Write(hist_name.c_str(), TObject::kOverwrite);

        f->cd();
        f->Write("", TObject::kOverwrite);
        std::cout << "[saved] " << rootfile << "  dir=" << hist_dir
                  << "  obj=" << hist_name << std::endl;
    }
}
