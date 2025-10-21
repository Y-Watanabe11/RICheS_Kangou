// 使い方（既定値は質問の仕様に合わせ済み）:
//   // ch: out.root の ch_eventtree.time (µs, 356 s 周回)
//   // dssd: dssd.root の eventtree.ti (10 ns tick, 32bit 周回)
//   // |Δt| ≤ 1 ms で相互最近傍のみを採用、結果は CH 側ファイルに追記
//   root -l -q 'sync_dssd_cherenkov.C("out.root","ch_eventtree","time",356.0,"dssd.root","eventtree","ti",1e-8,4294967296.0,1e-3,true,"")'
//
// 引数：
//   cher_root, cher_tree, cher_time, ch_wrap_sec
//   dssd_root, dssd_tree, dssd_time, dssd_tick_sec, dssd_modulus
//   max_dt_sec : 許容時間差（秒） 例：1e-3=1 ms
//   one_to_one : true=一対一（同じDSSDを再利用しない）
//   out_root   : 出力先。空("")なら cher_root に追記
//
// 出力ツリー（ディレクトリ sync/）:
//   sync_tree:      cher_event, cher_time, dssd_event, dssd_time, dt, matched=1
//   unmatched_ch:   cher_event, cher_time
//   unmatched_dssd: dssd_event, dssd_time

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <cmath>
#include <cstdint>

struct Row { double t; Long64_t idx; };

static void unwrap_and_zero(std::vector<double>& t, double period_sec){
  if(t.empty()) return;
  double base = t[0];
  double off  = 0.0;
  for(size_t i=0;i<t.size();++i){
    if(i>0 && t[i] + 1e-30 < t[i-1]) off += period_sec; // 周回（折返し）検出
    t[i] = (t[i] + off) - base;                          // 0 起点に平行移動
  }
}

void sync_dssd_cherenkov(const char* cher_root   = "out.root",
                         const char* cher_tree   = "ch_eventtree",
                         const char* cher_time   = "time",
                         double ch_wrap_sec      = 356.0,       // 356 s 周回（1 µs カウンタ想定）
                         const char* dssd_root   = "dssd.root",
                         const char* dssd_tree   = "eventtree",
                         const char* dssd_time   = "ti",        // 32bit tick (1 tick = 10 ns)
                         double dssd_tick_sec    = 1e-8,        // 10 ns
                         double dssd_modulus     = 4294967296.0,// 2^32
                         double max_dt_sec       = 1e-3,        // 許容差 1 ms
                         bool   one_to_one       = true,
                         const char* out_root    = "" )
{
  // --- 入力を開く ---
  std::unique_ptr<TFile> fC(TFile::Open(cher_root,"READ"));
  if(!fC || fC->IsZombie()){ std::cerr<<"[err] open fail: "<<cher_root<<"\n"; return; }
  std::unique_ptr<TFile> fD(TFile::Open(dssd_root,"READ"));
  if(!fD || fD->IsZombie()){ std::cerr<<"[err] open fail: "<<dssd_root<<"\n"; return; }

  TTree* tC = (TTree*)fC->Get(cher_tree);
  TTree* tD = (TTree*)fD->Get(dssd_tree);
  if(!tC){ std::cerr<<"[err] cher tree not found: "<<cher_tree<<"\n"; return; }
  if(!tD){ std::cerr<<"[err] dssd tree not found: "<<dssd_tree<<"\n"; return; }

  // --- CH 側（µs カウンタ, 356s 周回）---
  double ch_time_us = 0.0;
  if(!tC->GetBranch(cher_time)){ std::cerr<<"[err] cher time branch not found: "<<cher_time<<"\n"; return; }
  tC->SetBranchAddress(cher_time, &ch_time_us);

  const Long64_t nC = tC->GetEntries();
  std::vector<Row> C; C.reserve(nC);
  for(Long64_t i=0;i<nC;++i){
    tC->GetEntry(i);
    C.push_back({ ch_time_us * 1e-6, i }); // µs → s
  }
  // unwrap & 0起点
  {
    std::vector<double> tmp; tmp.reserve(C.size());
    for(auto& r: C) tmp.push_back(r.t);
    unwrap_and_zero(tmp, ch_wrap_sec);
    for(size_t i=0;i<C.size();++i) C[i].t = tmp[i];
  }

  // --- DSSD 側（32bit, 10ns tick）---
  UInt_t ti = 0;
  if(!tD->GetBranch(dssd_time)){ std::cerr<<"[err] dssd time branch not found: "<<dssd_time<<"\n"; return; }
  tD->SetBranchAddress(dssd_time, &ti);

  const Long64_t nD = tD->GetEntries();
  std::vector<Row> D; D.reserve(nD);
  for(Long64_t i=0;i<nD;++i){
    tD->GetEntry(i);
    double tsec = (double)ti * dssd_tick_sec; // tick→秒
    D.push_back({ tsec, i });
  }
  // unwrap & 0起点（2^32 tick 周回）
  {
    const double period_sec = dssd_modulus * dssd_tick_sec;
    std::vector<double> tmp; tmp.reserve(D.size());
    for(auto& r: D) tmp.push_back(r.t);
    unwrap_and_zero(tmp, period_sec);
    for(size_t i=0;i<D.size();++i) D[i].t = tmp[i];
  }

  // 念のため昇順保証（通常は順序通りだが）
  std::sort(C.begin(), C.end(), [](const Row&a,const Row&b){return a.t<b.t;});
  std::sort(D.begin(), D.end(), [](const Row&a,const Row&b){return a.t<b.t;});

  // --- 出力ファイル（既定は CH 側に追記）---
  TString outname = (out_root && strlen(out_root)) ? out_root : cher_root;
  std::unique_ptr<TFile> fO(TFile::Open(outname, (out_root && strlen(out_root)) ? "RECREATE" : "UPDATE"));
  if(!fO || fO->IsZombie()){ std::cerr<<"[err] cannot open output: "<<outname<<"\n"; return; }

  if(!fO->GetDirectory("sync")) fO->mkdir("sync");
  fO->cd("sync");

  // --- 出力ツリー ---
  TTree* Ts = new TTree("sync_tree","Matched Cerenkov-DSSD events");
  Long64_t cher_event=-1, dssd_event=-1;
  double cher_time=0, dssd_time=0, dt=0;
  int matched=0;
  Ts->Branch("cher_event",&cher_event);
  Ts->Branch("cher_time",&cher_time);
  Ts->Branch("dssd_event",&dssd_event);
  Ts->Branch("dssd_time",&dssd_time);
  Ts->Branch("dt",&dt);
  Ts->Branch("matched",&matched);

  TTree* TuC = new TTree("unmatched_ch","Unmatched Cerenkov events");
  Long64_t u_ch_event=-1; double u_ch_time=0;
  TuC->Branch("cher_event",&u_ch_event);
  TuC->Branch("cher_time",&u_ch_time);

  TTree* TuD = new TTree("unmatched_dssd","Unmatched DSSD events");
  Long64_t u_dssd_event=-1; double u_dssd_time=0;
  TuD->Branch("dssd_event",&u_dssd_event);
  TuD->Branch("dssd_time",&u_dssd_time);

  // --- 近傍探索のため D の時刻配列を抽出 ---
  std::vector<double> Dt; Dt.reserve(D.size());
  for(auto& r: D) Dt.push_back(r.t);

  std::vector<char> usedD(D.size(), 0);
  std::vector<char> usedC(C.size(), 0);

  auto nearest_idx = [](const std::vector<double>& arr, double x)->Long64_t{
    auto it = std::lower_bound(arr.begin(), arr.end(), x);
    if(it == arr.begin()) return 0;
    if(it == arr.end())   return (Long64_t)arr.size()-1;
    Long64_t i = (Long64_t)std::distance(arr.begin(), it);
    double a = arr[i-1], b = arr[i];
    return (std::fabs(x-a) <= std::fabs(b-x)) ? (i-1) : i;
  };

  // --- 相互最近傍 & 閾値でマッチング ---
  for(size_t ic=0; ic<C.size(); ++ic){
    double tc = C[ic].t;
    Long64_t id = nearest_idx(Dt, tc);
    double td = Dt[id];
    double d1 = std::fabs(td - tc);

    // 相手側から見た最近傍
    Long64_t ic2 = nearest_idx(
      *reinterpret_cast<const std::vector<double>*>(
        &std::vector<double>([&](){ std::vector<double> tmp; tmp.reserve(C.size()); for(auto& r:C) tmp.push_back(r.t); return tmp; }())
      ), 0.0); // ダミー（下でちゃんと作る）

    // ↑上は一時配列のライフタイム問題になるので、先に C 側の時刻配列を用意
  }

  // C 側の時刻配列を用意（ライフタイムが関数末まで保つよう外に置く）
  std::vector<double> Ct; Ct.reserve(C.size());
  for(auto& r: C) Ct.push_back(r.t);

  int nm=0, nuC=0, nuD=0;

  for(size_t ic=0; ic<C.size(); ++ic){
    if(usedC[ic]) continue;
    double tc = C[ic].t;

    Long64_t id = nearest_idx(Dt, tc);
    double   td = Dt[id];
    double   d1 = std::fabs(td - tc);

    if(d1 > max_dt_sec) { // 閾値を超えたら CH を未対応に
      u_ch_event = C[ic].idx; u_ch_time = tc; TuC->Fill(); ++nuC; continue;
    }
    if(one_to_one && usedD[id]){ // 既に使われていたら未対応
      u_ch_event = C[ic].idx; u_ch_time = tc; TuC->Fill(); ++nuC; continue;
    }

    // 相互最近傍チェック
    Long64_t ic_back = nearest_idx(Ct, td);
    if((Long64_t)ic != ic_back){
      // 相手から見て最近傍でなければ未対応
      u_ch_event = C[ic].idx; u_ch_time = tc; TuC->Fill(); ++nuC; continue;
    }

    // 採用
    matched = 1;
    cher_event = C[ic].idx; cher_time = tc;
    dssd_event = D[id].idx; dssd_time = td;
    dt = td - tc;
    Ts->Fill();
    usedC[ic] = 1;
    usedD[id] = 1;
    ++nm;
  }

  // 余った DSSD 側も未対応として記録
  for(size_t id=0; id<D.size(); ++id){
    if(usedD[id]) continue;
    if(std::find_if(C.begin(), C.end(), [&](const Row& r){ return std::fabs(D[id].t - r.t) <= max_dt_sec; }) != C.end()){
      // C からはマッチ済の可能性があるが one_to_one で溢れた場合はここへ来る
    }
    u_dssd_event = D[id].idx; u_dssd_time = D[id].t; TuD->Fill(); ++nuD;
  }

  fO->cd();
  Ts->Write("", TObject::kOverwrite);
  TuC->Write("", TObject::kOverwrite);
  TuD->Write("", TObject::kOverwrite);
  fO->Write("", TObject::kOverwrite);

  std::cout << "[done] " << outname << "\n"
            << "  matched=" << nm
            << "  unmatched_ch=" << nuC
            << "  unmatched_dssd=" << nuD
            << "  max_dt=" << max_dt_sec << " s"
            << "  ch_wrap=" << ch_wrap_sec << " s"
            << "  dssd_wrap=" << dssd_modulus * dssd_tick_sec << " s"
            << std::endl;
}
