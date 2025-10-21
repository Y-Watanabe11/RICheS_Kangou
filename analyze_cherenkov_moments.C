// analyze_cherenkov_moments.C
// 使い方（例）:
//   root -l -q 'analyze_cherenkov_moments.C("out.root", /*min_photons=*/50, /*flipY=*/true, /*save_canv=*/true, "ring_canv", /*pix_pitch_mm=*/6.0)'
//
// 引数：
//   rootfile      : 入出力ROOT（UPDATE）
//   min_photons   : 総光子数の下限（暗電流除去）
//   flipY         : y=8 を上にするか（既存ヒストと整合）
//   save_canv     : 診断キャンバス保存フラグ（ring_canv/）
//   canv_dir      : キャンバス保存ディレクトリ
//   pix_pitch_mm  : ピクセルピッチ(mm)。<=0 なら mm 変換をしない
//
// 出力：
//   TTree "ch_ringtree"（上書き）に、event, sumQ, max_ch, max_q,
//   max_ix, max_iy（flip適用後）, xc, yc, rc（px）, xc_mm, yc_mm（mm）,
//   a, b（px; 主/副軸スケール）, ok（=1:計算成功）
//
// 変更点：
//   - 画像モーメント（M00,M10,M01,u20,u02,u11）は 3×3 ではなく **全8×8画素**で計算。
//   - flipY は保存・描画の座標系にのみ影響（元の意図を踏襲）。

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TString.h>
#include <TROOT.h>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <cctype>

static TDirectory* ensure_dir(TFile* f, const char* name){
    if(!f) return nullptr;
    TDirectory* d = (TDirectory*)f->Get(name);
    if(!d) d = f->mkdir(name);
    return d;
}
static inline void ch_to_ixiy(int ch, int& ix, int& iy){
    int idx = ch-1;
    ix = (idx % 8) + 1;  // 1..8
    iy = (idx / 8) + 1;  // 1..8
}

void analyze_cherenkov_moments(const char* rootfile="out.root",
                               double min_photons=50.0,
                               bool flipY=true,
                               bool save_canv=true,
                               const char* canv_dir="ring_canv",
                               double pix_pitch_mm=6.0)
{
    gROOT->SetBatch(kTRUE);

    std::unique_ptr<TFile> f(TFile::Open(rootfile,"UPDATE"));
    if(!f || f->IsZombie()){ std::cerr << "[err] open fail: " << rootfile << "\n"; return; }

    TTree* t = (TTree*)f->Get("ch_eventtree");
    if(!t){ std::cerr << "[err] ch_eventtree not found\n"; return; }

    // ブランチ接続
    double ch[64] = {0};
    for(int i=0;i<64;++i){
        if(!t->GetBranch(Form("ch%d", i+1))){
            std::cerr<<"[err] missing branch ch"<<(i+1)<<"\n"; return;
        }
        t->SetBranchAddress(Form("ch%d", i+1), &ch[i]);
    }

    // 出力ツリー（上書き）
    if (f->Get("ch_ringtree")) f->Delete("ch_ringtree;*");
    TTree* rt = new TTree("ch_ringtree","Cherenkov ring by ALL-pixel image-moment");

    Long64_t event=0;
    double sumQ=0;
    int max_ch=0, max_ix=0, max_iy=0;
    double max_q=0;

    double xc=0, yc=0, rc=0;       // px（重心 & 半径）
    double xc_mm=0, yc_mm=0;       // mm
    double a=0, b=0;               // 主/副軸スケール（px）
    int ok=0;

    rt->Branch("event",&event);
    rt->Branch("sumQ",&sumQ);
    rt->Branch("max_ch",&max_ch);
    rt->Branch("max_q",&max_q);
    rt->Branch("max_ix",&max_ix);
    rt->Branch("max_iy",&max_iy);
    rt->Branch("xc",&xc);
    rt->Branch("yc",&yc);
    rt->Branch("rc",&rc);
    rt->Branch("a",&a);
    rt->Branch("b",&b);
    rt->Branch("xc_mm",&xc_mm);
    rt->Branch("yc_mm",&yc_mm);
    rt->Branch("ok",&ok);

    TDirectory* dcanv = save_canv ? ensure_dir(f.get(), canv_dir) : nullptr;

    const int NX=8, NY=8;
    const Long64_t nent = t->GetEntries();
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(50);
    gStyle->SetPaintTextFormat("g");

    for(event=0; event<nent; ++event){
        t->GetEntry(event);

        // 8x8 配列＆合計・最大
        double img[NY][NX];
        sumQ = 0; max_q = -1e300; max_ch = 1;
        for(int iy=0; iy<NY; ++iy){
            for(int ix=0; ix<NX; ++ix){
                double v = ch[iy*NX + ix];
                img[iy][ix] = v;
                sumQ += v;
                if(v > max_q){ max_q = v; max_ch = iy*NX + ix + 1; }
            }
        }
        if(sumQ < min_photons){
            // 暗電流とみなしてスキップ（Fillはせずに次へ）
            continue;
        }

        // 最明画素（1..8）→ flip 系に直すのは最後（メタ情報用）
        int p_ix, p_iy;
        ch_to_ixiy(max_ch, p_ix, p_iy); // 1..8（下から上が 1..8）

        // ---- 画像モーメント（全 8×8 画素）----
        // ゼロ次・一次モーメント
        double M00=0, M10=0, M01=0;
        for(int iy=1; iy<=NY; ++iy){
            for(int ix=1; ix<=NX; ++ix){
                double w = img[iy-1][ix-1];
                if(w <= 0) continue;
                double yy = flipY ? (NY - iy + 1) : iy; // 保存系のy
                M00 += w;
                M10 += w * ix;
                M01 += w * yy;
            }
        }

        ok = 0; xc=yc=rc=a=b=xc_mm=yc_mm=0;
        if(M00 > 0){
            // 重心
            xc = M10 / M00;
            yc = M01 / M00;

            // 二次中心モーメント（全 8×8）
            double u20=0, u02=0, u11=0;
            for(int iy=1; iy<=NY; ++iy){
                for(int ix=1; ix<=NX; ++ix){
                    double w = img[iy-1][ix-1];
                    if(w <= 0) continue;
                    double yy = flipY ? (NY - iy + 1) : iy;
                    double dx0 = ix - xc;
                    double dy0 = yy - yc;
                    u20 += w * dx0*dx0;
                    u02 += w * dy0*dy0;
                    u11 += w * dx0*dy0;
                }
            }

            // 共分散（正規化）
            double cxx = u20 / M00;
            double cyy = u02 / M00;
            double cxy = u11 / M00;

            // 固有値（λ1 ≥ λ2）
            double tr  = cxx + cyy;
            double det = cxx*cyy - cxy*cxy;
            double disc = tr*tr - 4*det;
            if(disc < 0) disc = 0;
            double s = std::sqrt(disc);
            double l1 = 0.5*(tr + s);
            double l2 = 0.5*(tr - s);
            if(l1 < l2) std::swap(l1,l2);

            // 軸スケール（1σ）
            a = (l1>0) ? std::sqrt(l1) : 0;
            b = (l2>0) ? std::sqrt(l2) : 0;

            // リング半径推定（主/副軸スケールの平均）
            rc = 0.5*(a + b);

            // 最明画素（保存系へ）
            max_ix = p_ix;
            max_iy = flipY ? (NY - p_iy + 1) : p_iy;

            // mm 変換（必要なら）
            if(pix_pitch_mm > 0){
                xc_mm = (xc - 0.5) * pix_pitch_mm;
                yc_mm = (yc - 0.5) * pix_pitch_mm;
            }

            ok = 1;
        }

        rt->Fill();

        // 診断キャンバス
        if(save_canv && dcanv && ok){
            dcanv->cd();
            TString cname = Form("evt_%06lld", event);
            if (dcanv->Get(cname)) dcanv->Delete(TString(cname+";*"));
            TCanvas* c = new TCanvas(cname, Form("Event %lld (sumQ=%.1f)", event, sumQ), 900, 820);

            TH2D h("h",";X bin;Y bin", NX, 0.5, NX+0.5, NY, 0.5, NY+0.5);
            for(int iy=0; iy<NY; ++iy){
                for(int ix=0; ix<NX; ++ix){
                    int ybin = flipY ? (NY - iy) : (iy + 1);
                    h.SetBinContent(ix+1, ybin, img[iy][ix]);
                }
            }
            h.SetMarkerSize(1.2);
            h.Draw("COLZ TEXT0");

            // モーメント楕円（1σ）
            TEllipse ell(xc, yc, a, b, 0, 360, 0);
            ell.SetLineColor(kRed+1);
            ell.SetLineWidth(3);
            ell.SetFillStyle(0);
            ell.Draw("same");

            // 中心マーカー
            TEllipse ctr(xc, yc, 0.12, 0.12); ctr.SetFillColor(kBlack); ctr.Draw("same");

            c->Update();
            c->Write(cname, TObject::kOverwrite);
            delete c;
        }
    }

    f->cd();
    rt->Write("", TObject::kOverwrite);
    f->Write("", TObject::kOverwrite);

    std::cout << "[done] wrote ch_ringtree"
              << (save_canv ? " and ring_canv/*" : "")
              << " (ALL-pixel image-moment method)"
              << "  file=" << rootfile
              << "  min_photons=" << min_photons
              << std::endl;
}
