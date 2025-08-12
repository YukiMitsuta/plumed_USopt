アンブレラサンプリングパラメータ最適化手法の使い方。

1．インストール
このコードはplumedに追加することで機能します。まずはサイトからplumedをダウンロードしてください。
現在のコードはplumed-2.9.0で機能することを確認しています。2．8以前のバージョンでは機能しません。
ダウンロードしたらsrc/biasのディレクトリにこのプログラムをコピーします。

$ cp plumed_USopt/src/bias/* /path/to/plumed-2.9.0/src/bias

また、このプログラムはｃ＋＋の行列計算ライブラリであるEigenを使用しているため、Eigenを事前にインストールしてください。
Ubuntuであればaptを使用してダウンロードできます。

$ sudo apt install libeigen3-dev

加えて、plumedをインストールするときにCPPFLAGSでEigenのライブラリを指定する必要があります。

$ ./configure --prefix=/installing/path/to/plumed-2.9.0 \
CPPFLAGS="-I/path/to/eigen3/include/eigen3/" \
（そのほかの設定）

＊Ubuntuのaptでインストールした場合のパスは/usr/include/eigen3/です。

2．チュートリアル（alanine dipeptide）
このコードはplumedでアンブレラサンプリングができる系であれば、どのプログラムでも利用可能です。
今回はチュートリアルのファイルとして、Gromacsでのalanine dipeptide計算を添付しました。

actionにはOPTIMIZERESTRAINTHOTELLINGを使用します。そのオプションは以下の通りです。
OFFDIAGONAL        バイアスポテンシャルの非対角項を入れるかどうか
ARG                バイアスをかけるCV
AT　　　　　　　  　　バイアスのreference pointの初期値
TARGET　　　　  　　 最適化のtarget point （the mean pointがこの位置になる用最適化される）
KAPPA　　　　　　  　バイアスポテンシャルの初期値。行列の次元に合わせてKAPPA0, KAPPA1...と指定する。
OPTSTRIDE          最適化のstride
CVSTEPSIZE         CVに対する最適化ステップサイズ（学習率）
KAPPASTEPSIZE      バイアスポテンシャルに対する最適化ステップサイズ（学習率）
OPTMETHOD          確率的最適化の手法
BETA1              確率的最適化のパラメータ（beta_1） 
BETA2              確率的最適化のパラメータ（beta_2）
EPSILON            確率的最適化のパラメータ（epsilon）
EXP_DECAYING_AVER  exponential decaying averageの項
HOTELLINGMAX　　　　Hotelling距離での上限（指定しても最適化が改善しないことが分かったので、99999をいれて無効化してください）↲
TARGETSIGMA　　　　 最適化の分散の上限
IGNORECV　　　　　　 最適化を無視する反応座標（無視する次元に1.0を入れてください）

TARGETSIGMAの次元はCVの次元ではなく分散共分散行列を固有値展開したときの固有値の大きな順番に対応することに注意が必要です。
IGNORECVは、CVに上限や下限があるとき、最適化が収束しなくなる場合があるので、その周辺では最適化を無視するよう設定してください（例えばcoordination numberなど）。

出力オプションは以下の通りです。label.arg_cntrもしくはlabel.arg1_arg2_kappaといった指定が必要であることに注意してください。
_cntr   最適化された参照点
_kappa  最適化されたバイアスポテンシャルのパラメータ
_mean   サンプリングの平均点
_Kgrad  バイアスポテンシャルの最適化勾配
_CVgrad 反応座標の最適化勾配

Step1:最適化計算を行う。
MD計算をする際の実行ファイルとしてcalljob.shを用意しました。必要に応じて修正してください。
最適化計算をする際のスクリプトは以下の部分です。

$ gmx grompp -f ./npt.mdp -o ./npt.tpr -c ./min.gro -p ./topol.top
$ gmx mdrun -deffnm npt -plumed ./plumed_npt.dat↲

Step2:最適化の結果を解析し、plumed.datを作成する。
計算結果で_cntrと_kappaをPRINTで出力させ、その最後の時間の値を解析すれば、最適化されたウィンドウのバイアスポテンシャルを指定できます。
ここでは、スクリプトのmkplumed.pyを用意したので、それを利用してください。
また、計算系を変えるときにはこのスクリプトの内容を適切に変更してください。

Step3:ウィンドウのサンプリングを行う。
plumed.datが作成できたら、後はメインの計算を実行すればアンブレラサンプリング計算ができます。

$ gmx grompp -f ./run.mdp -o ./run.tpr -c ./npt.gro -p ./topol.top
$ gmx mdrun -deffnm run -plumed ./plumed.dat
