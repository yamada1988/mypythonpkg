# mypythonpkg
python packages for HPC

# Openmm
myopenmm/MDCond.py
.inpファイルを読み込んでmdrunを行うメソッドが書かれたスクリプト。

incropenmm/xtcreporter.py
結果をxtcファイル形式で書き出すスクリプト。

incropenmm/mdrun.py
実際にmdrunを行うスクリプト。

使い方：

incropenmm/ディレクトリをクローン

PYTHONPATHにincropenmm/ディレクトリを追加

.inpファイルを作成、各種パラメータ設定を行う

incropenmm/mdrun.pyが使用例スクリプト。
