⓪ define_path.m
pathを自分のものに設定する

① save_dtacqdata.m
指定実験日の磁気プローブデータをmat形式で
pathname.rawdata/date/に保存(一度やればok)

② main_pcb200ch_auto.m
実験ログからdtacqショット番号などを「自動取得」して
①で保存しておいたデータによる磁気面をプロット

③ main_pcb200ch_manual.m
dtacqショット番号などを「手動で設定」して
①で保存しておいたデータによる磁気面をプロット

④ main_ionflow_pcb200ch_auto.m
実験ログからICCD撮影パラメータ、dtacqショット番号などを「自動取得」して
イオンフロー、①で保存しておいたデータによる磁気面をプロット

⑤ main_ionflow_manual.m
ICCD撮影パラメータを「手動で設定」して
イオンフローをプロット

⑥ main_ionvdist_pcb200ch_auto.m
実験ログからICCD撮影パラメータ、dtacqショット番号などを「自動取得」して
イオン速度分布、イオンフロー、①で保存しておいたデータによる磁気面をプロット