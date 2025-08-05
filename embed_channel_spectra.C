// 使い方：
//   root -l -q 'embed_channel_spectra.C()'
//   root -l -q 'embed_channel_spectra.C("out.root","hist_spectrum",true,1,64,512,0,4095,false,0)'

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TDirectory.h>
#include <TString.h>
#include <TROOT.h>
#include <TExec.h>          // ← 追加
#include <iostream>
#include <vector>
#include <limits>
#include <memory>
#include <cctype>
#include <cmath>

static TDirectory* ensure_dir(TFile* f, const char* dname){
    if (!f) return nullptr;
    TDirectory* d = (TDirectory*)f->Get(dname);
    if (!d) d = f->mkdir(dname);
    return d;
}
static std::string trim(const std::string& s){
    auto a=s.begin(), b=s.end();
    while (a!=b && std::isspace((unsigned char)*a)) ++a;
    while (b!=a && std::isspace((unsigned char)*(b-1))) --b;
    std::string t(a,b);
    if (!t.empty() && (t.front()=='"' || t.front()=='\'')) t.erase(t.begin());
    if (!t.empty() && (t.back()=='"'  || t.back()=='\'' )) t.pop_back();
    return t;
}
static inline bool isnum(double x){ return std::isfinite(x); }

void embed_channel_spectra(const char* rootfile_cstr="",
                           const char* out_dir="hist_spectrum",
                           bool overwrite=true,
                           int ch_first=1, int ch_last=64,
                           int nbins=256,
                           double xmin=NAN, double xmax=NAN,
                           bool autoscan=true,
                           int compr_level=-1)
{
    // 対話でROOTファイル取得
    std::string rootfile = trim(rootfile_cstr ? std::string(rootfile_cstr) : std::string());
    if (rootfile.empty()){
        std::cout << "ROOTファイルのパスを入力してください（ドラッグ&ドロップ可）: ";
        std::getline(std::cin, rootfile);
        rootfile = trim(rootfile);
        if (rootfile.empty()){ std::cerr << "[err] ROOTファイル未指定\n"; return; }
    }

    // バッチモード（描画抑制）
    gROOT->SetBatch(kTRUE);

    std::unique_ptr<TFile> f(TFile::Open(rootfile.c_str(), "UPDATE"));
    if (!f || f->IsZombie()){ std::cerr << "[err] open fail: " << rootfile << "\n"; return; }

    int old_compr = f->GetCompressionLevel();
    if (compr_level >= 0) f->SetCompressionLevel(compr_level);

    // 入力ツリー
    TTree* t = (TTree*)f->Get("ch_eventtree");
    if (!t){ std::cerr << "[err] ch_eventtree not found\n"; return; }

    // 対象chチェック
    if (ch_first < 1) ch_first = 1;
    if (ch_last  > 64) ch_last = 64;
    if (ch_first > ch_last){ std::cerr << "[err] invalid channel range\n"; return; }

    // ブランチ接続
    double ch[64] = {0};
    for (int i=0;i<64;++i){
        auto* br = t->GetBranch(Form("ch%d", i+1));
        if (!br){ std::cerr << "[err] missing branch: ch" << (i+1) << "\n"; return; }
        t->SetBranchAddress(Form("ch%d", i+1), &ch[i]);
    }

    const Long64_t nent = t->GetEntries();
    if (nent <= 0){ std::cerr << "[err] empty tree\n"; return; }

    // 出力ディレクトリ
    TDirectory* dd = ensure_dir(f.get(), out_dir);
    if (!dd){ std::cerr << "[err] cannot get/make dir: " << out_dir << "\n"; return; }

    // 範囲決定
    bool use_global_range = isnum(xmin) && isnum(xmax) && (xmin < xmax) && !autoscan;
    std::vector<double> vmin(65, +std::numeric_limits<double>::infinity());
    std::vector<double> vmax(65, -std::numeric_limits<double>::infinity());

    if (!use_global_range){
        for (Long64_t ev=0; ev<nent; ++ev){
            t->GetEntry(ev);
            for (int c=ch_first; c<=ch_last; ++c){
                double v = ch[c-1];
                if (!isnum(v)) continue;
                if (v < vmin[c]) vmin[c] = v;
                if (v > vmax[c]) vmax[c] = v;
            }
            if ((ev % 200000)==0 && ev>0) std::cout << "[scan] " << ev << "/" << nent << "\n";
        }
        for (int c=ch_first; c<=ch_last; ++c){
            if (!std::isfinite(vmin[c]) || !std::isfinite(vmax[c]) || vmin[c]==vmax[c]){
                vmin[c] = 0.0; vmax[c] = 1.0;
            } else {
                double span = vmax[c]-vmin[c];
                vmin[c] -= 0.01*span; vmax[c] += 0.01*span;
            }
        }
    }

    // ヒスト作成（chごと）
    std::vector<TH1D*> hs(65, nullptr);
    for (int c=ch_first; c<=ch_last; ++c){
        TString hname = Form("spec_ch%02d", c);
        if (!overwrite && dd->Get(hname)) {
            std::cout << "[skip] exist: " << hname << "\n";
            continue;
        }
        double lo = use_global_range ? xmin : vmin[c];
        double hi = use_global_range ? xmax : vmax[c];
        if (!(hi>lo)) { lo = 0.0; hi = 1.0; }

        TH1D* h = new TH1D(hname, Form("Spectrum of ch%d;Value;Counts", c), nbins, lo, hi);
        h->SetDirectory(nullptr);
        h->Sumw2(false);

        // --- ここからログ表示の埋め込み ---
        h->SetMinimum(0.5);  // ゼロビンでもlog表示可能に（適宜変更可）
        TExec* ex = new TExec(Form("logy_ch%02d", c), "if(gPad){ gPad->SetLogy(1); }");
        h->GetListOfFunctions()->Add(ex);   // 描画時に自動でlogYへ
        // --- ここまで ---

        hs[c] = h;
    }

    // フィル
    for (Long64_t ev=0; ev<nent; ++ev){
        t->GetEntry(ev);
        for (int c=ch_first; c<=ch_last; ++c){
            TH1D* h = hs[c];
            if (!h) continue;
            double v = ch[c-1];
            if (!isnum(v)) continue;
            h->Fill(v);
        }
        if ((ev % 200000)==0 && ev>0) std::cout << "[fill] " << ev << "/" << nent << "\n";
    }

    // 書き込み
    int nwritten = 0, nskip = 0;
    for (int c=ch_first; c<=ch_last; ++c){
        TH1D* h = hs[c];
        if (!h) { ++nskip; continue; }
        dd->WriteTObject(h, h->GetName(), "Overwrite");
        delete h; hs[c] = nullptr;
        ++nwritten;
    }

    f->Write("", TObject::kOverwrite);
    if (compr_level >= 0) f->SetCompressionLevel(old_compr);

    std::cout << "[done] file=" << rootfile
              << "  channels=" << (ch_last - ch_first + 1)
              << "  wrote=" << nwritten
              << "  skipped=" << nskip
              << "  entries=" << nent
              << "  out_dir=" << out_dir
              << (use_global_range ? "  range=global" : "  range=autoscan(per-channel)")
              << "  (Y-axis log enabled)"
              << std::endl;
}
