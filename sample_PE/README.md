Polyethylene(PE)単成分系でincrOpenMM計算を行うインプットファイル群が置かれている。

inpdir/
MD計算のためのインプットファイルである.inpファイル置き場。

PE_0XX/ 先頭の溶質PEに関してXXmer部分までCore, XXmer以降がGhost原子の系のMD計算作業ディレクトリ。
MD/ トラジェクトリ(.xtc)およびログ(.log)が書き出される。
OUT/ ジョブスクリプトの.outおよび.errが書き出される。
SYS/ 系のトポロジー（.top, .itp）および構造ファイル(.gro)置き場。system*.groはMD計算後に更新される。
script/ MD計算を行うsample.pyおよびエネルギー最小化計算を行うmin.pyが置かれている。
system.ndx 溶質ポリマーのCore, Ghost部分を指定するファイル。
mdrun_02.sh ../inpdir/stage2/npt_ghost.inpをインプットとしてMD計算を行うスクリプト。
min.sh  ../inpdir/stage0/min.inpをインプットとしてエネルギー最小化計算を行うスクリプト。
min.sh => mdrun_02.shの順でテストMD計算を行うことができる。
