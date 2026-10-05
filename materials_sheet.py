"""
Export the consolidated material-property database from sims.py.

This script creates a single CSV file containing every workpiece, tool, and
base material with a common set of columns. You can then import the CSV into
Google Sheets (File > Import > Upload > Replace current sheet).

Optional automated upload:
    1. Create a Google Cloud service account and download its JSON key.
    2. Share your target Google Sheet with the service account email.
    3. Set GOOGLE_SHEETS_CREDENTIALS to the path of the JSON key.
    4. python -m pip install gspread google-auth
    5. python materials_sheet.py <SPREADSHEET_ID> <WORKSHEET_NAME>

Usage:
    python materials_sheet.py
    python materials_sheet.py 1SgXsug8FHFSZo2Z1j3qCiwfU8V8IxjhMgHmVqg2l34k Sheet1
"""

import csv
import os
import sys

from sims import build_materials_export_table


LOCAL_CSV = "materials_export.csv"


def to_csv(rows, filepath):
    """Write a list of dicts to CSV, using the dict keys as the header."""
    if not rows:
        print("No material data to export.")
        return

    fieldnames = list(rows[0].keys())
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved {len(rows)} materials to {filepath}")
    print("Columns:", ", ".join(fieldnames))


def upload_to_google_sheets(rows, spreadsheet_id, worksheet_name):
    """Upload rows to a Google Sheet using gspread with a service account."""
    try:
        import gspread
        from google.oauth2.service_account import Credentials
    except ImportError:
        print(
            "\nAutomated upload needs gspread. Install it with:\n"
            "    python -m pip install gspread google-auth\n"
            "Then run this script again with a spreadsheet ID and worksheet name."
        )
        return False

    creds_path = os.environ.get("GOOGLE_SHEETS_CREDENTIALS")
    if not creds_path:
        print(
            "\nSet GOOGLE_SHEETS_CREDENTIALS to the path of your service-account JSON key.\n"
            "Example:\n"
            "    export GOOGLE_SHEETS_CREDENTIALS=/path/to/service-account.json"
        )
        return False

    scopes = [
        "https://www.googleapis.com/auth/spreadsheets",
        "https://www.googleapis.com/auth/drive",
    ]
    creds = Credentials.from_service_account_file(creds_path, scopes=scopes)
    client = gspread.authorize(creds)

    sheet = client.open_by_key(spreadsheet_id)
    try:
        worksheet = sheet.worksheet(worksheet_name)
    except gspread.WorksheetNotFound:
        worksheet = sheet.add_worksheet(title=worksheet_name, rows="100", cols="20")

    # Clear existing content and write header + data
    worksheet.clear()
    fieldnames = list(rows[0].keys())
    values = [fieldnames] + [[str(row.get(col, "")) for col in fieldnames] for row in rows]
    worksheet.update(values)

    print(f"\nUploaded {len(rows)} materials to Google Sheet {spreadsheet_id}, worksheet '{worksheet_name}'.")
    return True


def main():
    rows = build_materials_export_table()

    # Always write a local CSV for easy manual import into Google Sheets.
    to_csv(rows, LOCAL_CSV)

    # Optional: upload directly if a spreadsheet ID was provided.
    if len(sys.argv) >= 2:
        spreadsheet_id = sys.argv[1]
        worksheet_name = sys.argv[2] if len(sys.argv) >= 3 else "Sheet1"
        upload_to_google_sheets(rows, spreadsheet_id, worksheet_name)
    else:
        print(
            "\nTo upload automatically, run:\n"
            "    python materials_sheet.py <SPREADSHEET_ID> [WORKSHEET_NAME]\n"
            "and set GOOGLE_SHEETS_CREDENTIALS to your service-account JSON key.\n"
            "Or manually import the CSV into Google Sheets: File > Import > Upload."
        )


if __name__ == "__main__":
    main()
