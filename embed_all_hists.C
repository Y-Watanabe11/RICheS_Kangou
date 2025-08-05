// 高速版: TH2Dを再利用・WriteTObjectで上書き・圧縮任意・コンパイル推奨
// 実行例（コンパイル実行が推奨）:
//   root -l -q 'embed_all_hists.C+("out.root","hist",true,true,0,-1,1,-1,-1,"g",1.4,false,"canv", -1)'
//
// 主要引数：
//   rootfile_cstr : 入出力ROOT（UPDATE、空だと対話入力）
//   hist_dir      : ヒスト保存先ディレクトリ
//   overwrite     : 上書き(true) / 既存維持(false)
//   flipY         : 画像的に上をy=8
//   first/last    : 対象event範囲（-1で末尾）
//   every         : 間引き（1=全件）
//   zmin/zmax     : z固定（<0で自動）
//   textfmt       : "g", ".0f", "4.1f" など（%は付けない）
//   textsize      : 数値文字サイズ（MarkerSize）
//   make_canv     : 数値付きCanvasも保存するか（既定false=しない）
//   canv_dir      : Canvas保存先
//   compr_level   : 圧縮レベル（-1=現状維持、0=無圧縮, 1〜9=圧縮）

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TExec.h>
#include <TROOT.h>
#include <iostream>
#include <memory>
#include <string>
#include <cctype>

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

void embed_all_hists(const char* rootfile_cstr="",
                     const char* hist_dir="hist",
                     bool overwrite=true,
                     bool flipY=true,
                     Long64_t first=0,
                     Long64_t last=-1,
                     int every=1,
                     double zmin=-1, double zmax=-1,
                     const char* textfmt="g",
                     double textsize=1.4,
                     bool make_canv=false,
                     const char* canv_dir="canv",
                     int compr_level=-1)
{
    // 入力
    std::string rootfile = trim(rootfile_cstr ? std::string(rootfile_cstr) : std::string());
    if (rootfile.empty()){
        std::cout << "埋め込み先の ROOT ファイルのパスを入力してください（ドラッグ&ドロップ可）: ";
        std::getline(std::cin, rootfile);
        rootfile = trim(rootfile);
        if (rootfile.empty()){ std::cerr << "[err] ROOTファイル未指定\n"; return; }
    }

    // バッチモード（描画抑制）
    gROOT->SetBatch(kTRUE);

    std::unique_ptr<TFile> f(TFile::Open(rootfile.c_str(), "UPDATE"));
    if (!f || f->IsZombie()){ std::cerr << "[err] open fail: " << rootfile << "\n"; return; }

    // 圧縮設定（任意）
    int old_compr = f->GetCompressionLevel();
    if (compr_level >= 0) f->SetCompressionLevel(compr_level);

    // 入力ツリー
    TTree* t = (TTree*)f->Get("ch_eventtree");
    if (!t){ std::cerr << "[err] ch_eventtree not found\n"; return; }

    // ch1..ch64 を接続
    double ch[64] = {0};
    for (int i=0;i<64;++i){
        auto* br = t->GetBranch(Form("ch%d", i+1));
        if (!br){ std::cerr << "[err] missing branch: ch" << (i+1) << "\n"; return; }
        t->SetBranchAddress(Form("ch%d", i+1), &ch[i]);
    }

    const int nx=8, ny=8;
    Long64_t nall = t->GetEntries();
    if (nall <= 0){ std::cerr << "[err] empty ch_eventtree\n"; return; }
    if (last < 0 || last >= nall) last = nall - 1;
    if (first < 0) first = 0;
    if (first > last){ std::cerr << "[err] invalid range\n"; return; }

    TDirectory* dh = ensure_dir(f.get(), hist_dir);
    if (!dh){ std::cerr << "[err] cannot get/make dir: " << hist_dir << "\n"; return; }

    TDirectory* dc = nullptr;
    if (make_canv){
        dc = ensure_dir(f.get(), canv_dir);
        if (!dc){ std::cerr << "[err] cannot get/make dir: " << canv_dir << "\n"; return; }
    }

    // ここから高速化の肝：ヒストを一度だけ作り、毎回 Reset/SetName して再利用
    TH2D h("hwork", "work;X bin;Y bin", nx, 0.5, nx+0.5, ny, 0.5, ny+0.5);
    h.SetDirectory(nullptr);               // ディレクトリ非所属（軽い）
    h.SetMarkerSize(textsize);
    h.SetMarkerColor(kBlack);
    h.SetDrawOption("COLZ TEXT0");         // 保存後も数値が出る
    if (zmin>=0) h.SetMinimum(zmin);
    if (zmax>=0) h.SetMaximum(zmax);

    // 文字の描画形式（Canvas保存時や後の再描画用）
    gStyle->SetPaintTextFormat(textfmt);
    gStyle->SetOptStat(0);

    Long64_t nwritten = 0, nskip_exist = 0;

    for (Long64_t ev = first; ev <= last; ++ev){
        if (every > 1 && ((ev-first) % every)) continue;

        const TString hname = Form("histtree_%06lld", ev);
        if (!overwrite && dh->Get(hname)) { ++nskip_exist; continue; }

        // エントリ取得
        t->GetEntry(ev);

        // 再利用ヒストを埋める
        h.Reset("ICES");                               // 既存内容を速やかにクリア
        h.SetTitle(Form("Event %lld;X bin;Y bin", ev));
        h.SetName(hname);                              // 名前をその都度付与
        for (int iy=0; iy<ny; ++iy){
            const int row = iy*nx;
            const int ybin = flipY ? (ny - iy) : (iy + 1);
            for (int ix=0; ix<nx; ++ix){
                h.SetBinContent(ix+1, ybin, ch[row + ix]);
            }
        }

        // ディレクトリ切替不要：WriteTObjectで直接書き込み（上書き）
        dh->WriteTObject(&h, hname, "Overwrite");

        // 必要なら Canvas も保存（遅くなるので既定はオフ）
        if (make_canv){
            TString cname = Form("evt_%06lld", ev);
            TCanvas c(cname, h.GetTitle(), 1000, 850);
            // 再描画時にも反映されるよう TExec でStyle設定
            TExec ex("apply_textfmt", Form("gStyle->SetPaintTextFormat(\"%s\"); gStyle->SetOptStat(0);", textfmt));
            ex.Draw();
            h.Draw(h.GetDrawOption());     // "COLZ TEXT0"
            c.Modified(); c.Update();
            dc->WriteTObject(&c, cname, "Overwrite");
        }

        ++nwritten;

        // 大量書き込み時は定期的にフラッシュ
        if ((nwritten % 2000) == 0) {
            f->Flush();
            std::cout << "[info] written " << nwritten << " (event " << ev << ")\n";
        }
    }

    // 変更反映＆圧縮設定を戻す
    f->Write("", TObject::kOverwrite);
    if (compr_level >= 0) f->SetCompressionLevel(old_compr);

    std::cout << "[done] file=" << rootfile
              << "  wrote=" << nwritten
              << "  skipped_exist=" << nskip_exist
              << "  range=[" << first << "," << last << "]"
              << "  step=" << every
              << "  hist_dir=" << hist_dir
              << (make_canv ? TString::Format("  canv_dir=%s", canv_dir).Data() : "")
              << "  compression=" << (compr_level>=0?compr_level:old_compr)
              << std::endl;
}
