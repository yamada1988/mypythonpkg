# https://xlsxwriter.readthedocs.io/index.html
# conda install XlsxWriter   
import xlsxwriter
from xlsxwriter.utility import xl_rowcol_to_cell, xl_range_abs

nrow =  1
nline = 64

xlsxname = 'DAT/mus.xlsx'
wb = xlsxwriter.Workbook(xlsxname)

Nst = 50
Ned = 50
for k in range(Nst, Ned-1, -1):
    dirname = 'MELT_{0:03d}'.format(k)
    mus = [['' for i in range(nrow)] for j in range(nline)]
    errs = ['' for i in range(nline)] 
    print('dir:', k)
    for i in range(1, nline+1):
        for j in range(1, nrow+1):
            print('sln:', i, 'refs:', j)
            fname = dirname + '/ERmods/ermod_sln{0:02d}_refs{1:02d}/slvfe.tt'.format(i, j)
            with open(fname, 'rt') as f:
                lines = [line.strip() for line in f ]
            try:
                mu = lines[-3].split()[1]
                err = lines[-3].split()[2]
                mus[i-1][j-1] = float(mu)
                errs[i-1] = float(err)
            except:
                mu = ''   
                mus[i-1] = mu
                errs[i-1] = mu
        print(mus)
        print(errs)

    ws = wb.add_worksheet("mus_{0:03d}".format(k))

    for i in range(200):
        ws.set_row(i, 12.5)
        ws.set_column(i, i, 15.0)

    # header
    rfs_row = [ 'refs={0:02d}'.format(line) for line in range(1, nrow+1) ]
#    header_row = [''] + rfs_row + [ 'Average', 'Stdev', '95%error']
    header_row = [''] + rfs_row + ['Error']
    sln_col = [ 'soln={0:02d}'.format(line) for line in range(1, nline+1)]
    header_col = [''] + sln_col + [ 'Average', 'Stdev', '95%error']
    alphatable = {0:'A', 1:'B', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I', 9:'J', 10:'K'}

    for j in range(0, len(header_row)):
        ws.write_string(0, j, header_row[j])

    for i in range(0, len(header_col)):
        ws.write_string(i, 0, header_col[i])

    for i in range(1, nline+1):
        for j in range(1, nrow+1):
            try:
                ws.write_number(i, j, mus[i-1][j-1])
            except:
                ws.write_string(i, j, '')

    for i in range(1, nline+1):
        try:
            ws.write_number(i, nrow+1, errs[i-1])
        except:
            ws.write_string(i, nrow+1, '')


    # Set Formula
    for j in range(1, nrow+1):
        str_ = alphatable[j] + '2:' + alphatable[j] + '{0:d}'.format(nline+1)
        ws.write_formula(nline+1, j, '=AVERAGE(' + str_ +')')
        ws.write_formula(nline+2, j, '=STDEV(' +str_ + ')')
        ws.write_formula(nline+3, j, '='+alphatable[j]+'{0:d}/SQRT({1:d})'.format(nline+3, nline))

#    for i in range(1, nline+1):
#        str_ = 'B{0:d}:'.format(i+1) + alphatable[nrow] + '{0:d}'.format(i+1)
#        ws.write_formula(i, nrow+2, '=AVERAGE('+ str_ +')')
#        ws.write_formula(i, nrow+3, '=STDEV(' + str_ + ')') 
#        ws.write_formula(i, nrow+4, '='+ alphatable[nrow+2] +'{0:d}/SQRT({1:d})'.format(i+1,nrow))

    props = {
    "type": "3_color_scale",
    "max_color": "#ff0000", # red
    "mid_color": "#ffffff", # white
    "min_color": "#0000ff", # blue
    "max_type": "num",
    "mid_type": "num",
    "min_type": "num",
    "max_value": 2.0,
    "mid_value": -3.0,
    "min_value":  -8.0}
    ws.conditional_format(1, 1, nline, nrow, props)


wb.close()
