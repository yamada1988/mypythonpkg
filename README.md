# mypythonpkg
python packages for HPC

# incrOpenmm
incropenmm/MDCond.py
.inpファイルを読み込んでmdrunを行うメソッドが書かれたスクリプト。

incropenmm/xtcreporter.py
結果をxtcファイル形式で書き出すスクリプト。

incropenmm/mdrun.py
実際にmdrunを行うスクリプト。

sample_PE/
Polyethylene(PE)単成分系でincrOpenMM計算を行うインプットファイル群が置かれている。

使い方：

incropenmm/ディレクトリをクローン

PYTHONPATHにincropenmm/ディレクトリを追加

.inpファイルを作成、各種パラメータ設定を行う

incropenmm/mdrun.pyが使用例スクリプト。

# ermod_mod
incrOpenMMでrefs系のermod計算を行うための修正ermodのプログラム置き場。

# chain-increment-manual
OpenMMおよびincrOpenMMのマニュアル。

# sample_PE
鎖長N=100のPolyethylene(PE)単成分系の計算サンプルディレクトリ。

# sample_PEPVDF
鎖長N=100のPolyethylene(PE)とPolyvinylidenedifluoide(PVDF)二成分系の計算サンプル。
