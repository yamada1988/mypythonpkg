import xlsxwriter
from xlsxwriter.utility import xl_rowcol_to_cell, xl_range_abs

nraw =  5
nline = 52

xlsxname = 'DAT/mus.xlsx'
wb = xlsxwriter.Workbook(xlsxname)

Nst = 50
Ned = 50
for k in range(Nst, Ned-1, -1):
    dirname = 'MELT_{0:03d}'.format(k)
    mus = [['' for i in range(nraw)] for j in range(nline)]
    errs = [['' for i in range(nraw)] for j in range(nline)] 
    print('dir:', k)
    for i in range(1, nline+1):
        for j in range(1, nraw+1):
            print('sln:', i, 'refs:', j)
            fname = dirname + '/ERmods/ermod_sln{0:02d}_refs{1:02d}/slvfe.tt'.format(i, j)
            with open(fname, 'rt') as f:
                lines = [line.strip() for line in f ]
            try:
                mu = lines[-3].split()[1]
                err = lines[-3].split()[2]
                mus[i-1][j-1] = float(mu)
                errs[i-1][j-1] = float(err)
            except:
                mu = ''   
                mus[i-1] = mu
                errs[i-1] = mu
        print(mus)

    ws = wb.add_worksheet("mus_{0:03d}".format(k))
    # header
    header_row = ['', 'refs=1', 'refs=2', 'refs=3', 'refs=4', 'refs=5', 'Average', 'Stdev', '95%error']
    header_col = ['', 'soln=1', 'soln=2', 'soln=3', 'soln=4', 'soln=5', 'soln=6', 'soln=7', 'soln=8', 'soln=9', 'soln=10', 
                      'soln=11', 'soln=12', 'soln=13', 'soln=14', 'soln=15', 'soln=16', 'soln=17', 'soln=18', 'soln=19', 'soln=20', 
                      'soln=21', 'soln=22', 'soln=23', 'soln=24', 'soln=25', 'soln=26', 'soln=27', 'soln=28', 'soln=29', 'soln=30', 
                      'soln=31', 'soln=32', 'soln=33', 'soln=34', 'soln=35', 'soln=36', 'siln=37', 'soln=38', 'soln=39', 'soln=40',
                      'soln=41', 'soln=42', 'soln=43', 'soln=44', 'soln=45', 'soln=46', 'siln=47', 'soln=48', 'soln=49', 'soln=50', 
                      'soln=51', 'soln=52', 'Average', 'Stdev', '95%error']
    alphatable = {0:'A', 1:'B', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I', 9:'J', 10:'K'}

    for j in range(0, 8+1):
        ws.write_string(0, j, header_row[j])

    for i in range(0, len(header_col)):
        ws.write_string(i, 0, header_col[i])

    for i in range(1, nline+1):
        for j in range(1, nraw+1):
            try:
                ws.write_number(i, j, mus[i-1][j-1])
            except:
                ws.write_string(i, j, '')

    # Set Formula
    for j in range(1, nraw+1):
        str_ = alphatable[j] + '2:' + alphatable[j] + '{0:d}'.format(nline+1)
        ws.write_formula(nline+1, j, '=AVERAGE(' + str_ +')')
        ws.write_formula(nline+2, j, '=STDEV(' +str_ + ')')
        ws.write_formula(nline+3, j, '='+alphatable[j]+'55/SQRT({0:d})'.format(nline))

    for i in range(1, nline+1):
        ws.write_formula(i, nraw+1, '=AVERAGE(B{0:d}:F{0:d})'.format(i+1))
        ws.write_formula(i, nraw+2, '=STDEV(B{0:d}:F{0:d})'.format(i+1))                   
        ws.write_formula(i, nraw+3, '=H{0:d}/SQRT({1:d})'.format(i+1,nraw))

    props = {
    "type": "2_color_scale",
    "max_color": "#FF0000", # tomato
    "min_color": "#00FF00", # blue
    "max_type": "formula",
    "min_type": "formula",
    "max_value": -8.0,
    "min_value":  2.0}
    ws.conditional_format(1, 1, nline, nraw, props)


wb.close()
