import gspread
import pandas as pd

gc = gspread.oauth()

def getTS6log() -> pd.DataFrame:
    sh = gc.open_by_url('https://docs.google.com/spreadsheets/d/1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw/edit#gid=0')
    worksheet = sh.worksheet("log")
    expected_headers = ['number']  # Replace with your actual headers
    df = pd.DataFrame(worksheet.get_all_records(expected_headers=expected_headers))
    return df

def searchlog(df: pd.DataFrame, node: str, pat) -> pd.DataFrame:
    if node in df.columns:
        return df[df[node] == pat]
    else:
        raise KeyError(f"{node} は DataFrame の列に存在しません。")