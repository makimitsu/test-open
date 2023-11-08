%実験ログをwebからxlsx形式でダウンロード
ID = '1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
url_name = strcat('https://docs.google.com/spreadsheets/d/',ID,'/export?format=xlsx');
websave('exp_log.xlsx',url_name);

