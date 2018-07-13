#write_xlsx.py
import xlsxwriter
from xlsxwriter.utility import xl_rowcol_to_cell, xl_range_abs

mus = arr = [['' for i in range(10)] for j in range(10)]
#print(mus)

xlsxname = 'DAT/mus.xlsx'
wb = xlsxwriter.Workbook(xlsxname)

for k in range(50, 30, -1):
    dirname = 'MELT_{0:03d}'.format(k)
    mus = arr = [['' for i in range(10)] for j in range(10)]
    for i in range(1, 11):
        for j in range(1, 11):
            fname = dirname + '/ERmods/ermod_sln{0:02d}_ref{1:02d}/slvfe.tt'.format(i, j)
            with open(fname, 'rt') as f:
                lines = [line.strip() for line in f ]
            mu = lines[-3].split()[1]
            mus[i-1][j-1] = float(mu)

    print(mus)

    ws = wb.add_worksheet("mus_{0:03d}".format(k))
    # header
    header_row = ['', 'refs=1', 'refs=2', 'refs=3', 'refs=4', 'refs=5', 'refs=6', 'refs=7', 'refs=8', 'refs=9', 'refs=10', 'Average', 'Stdev', '95%error']
    header_col = ['', 'soln=1', 'soln=2', 'soln=3', 'soln=4', 'soln=5', 'soln=6', 'soln=7', 'soln=8', 'soln=9', 'soln=10', 'Average', 'Stdev', '95%error']
    alphatable = {0:'A', 1:'B', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I', 9:'J', 10:'K'}

    for j in range(0, 14):
        ws.write_string(0, j, header_row[j])

    for i in range(0, 14):
        ws.write_string(i, 0, header_col[i])

    for i in range(1, 11):
        for j in range(1, 11):
            ws.write_number(i, j, mus[i-1][j-1])

    for i in range(1, 11):
        ws.write_formula(i, 11, '=AVERAGE(B{0:d}:K{1:d})'.format(i+1, i+1))
        ws.write_formula(i, 12, '=STDEV(B{0:d}:K{1:d})'.format(i+1, i+1))

    for j in range(1, 11):
        A = alphatable[i]
        ws.write_formula(11, j, '=AVERAGE('+A+'2:'+A+'11)')
        ws.write_formula(12,  j, '=STDEV('+A+'2:'+A+'11)')

#    max_cell = xl_rowcol_to_cell(0, 7, row_abs=True, col_abs=True)
#    min_cell = xl_rowcol_to_cell(1, 7, row_abs=True, col_abs=True)
    props = {
    "type": "2_color_scale",
    "max_color": "#FF0000", # tomato
    "min_color": "#00FF00", # lightyellow
    "max_type": "formula",
    "min_type": "formula",
    "max_value": -6.0,
    "min_value": -1.0}
    ws.conditional_format(1, 1, 10, 10, props)

    ws.write_formula(12, 11, '=AVERAGE(B2:K11)')
    ws.write_formula(12, 12, '=STDEV(B2:K11)')
    ws.write_formula(12, 13, '=STDEV(B2:K11)*2.0/sqrt(100)')

wb.close()
