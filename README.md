# mypythonpkg
python packages for HPC

# Openmm
myopenmm/MDCond.py
.inpファイルを読み込んでmdrunを行うメソッドが書かれたスクリプト。

myopenmm/xtcreporter.py
結果をxtcファイル形式で書き出すスクリプト。

myopenmm/mdrun.py
実際にmdrunを行うスクリプト。

使い方：

myopenmm/ディレクトリをクローン

PYTHONPATHにmyopenmm/ディレクトリを追加

openmm.inpファイルを作成、各種パラメータ設定を行う

myopenmm/mdrun.pyが使用例スクリプト。
