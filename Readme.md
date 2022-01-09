User's Guide for fof

                       2022/01/09 牧野


# はじめに

この文章では、 FDPS を使った constant link での Friend-of-friend アル
ゴリズムによる粒子クラスタリングを行う fof について述べる。

# コンパイル方法等

基本的に FDPS の sample/c++/nbody の下にある nbody.cpp と同じだが、
PS_PATH を適切に (fdps の src があるディレクトリをさすように)設定す
る必要がある。

```
 make 
```
で実行ファイル fof が生成される。デフォルトでは OpenMP 並列版が生成さ
れる。Makefile をサンプルプログラム nbody.cpp の場合と同様に編集するこ
とで MPI 並列化も可能である。(テストしてないけど)
 

# 実行

実行時オプションについて説明する

 -i: 入力ファイルを指定する。入力ファイルは標準的な FDPS のアスキー形式で、

* ヘッダに時刻と粒子数がある
* 1行には、粒子の id番号、質量、座標の3成分がこの順である。

ものを想定している。それ以外のデータがあってもよいが読み飛ばす。つまり

```
time
n_particles
id[0]  mass[0] x[0] y[0] z[0] ...
id[1]  mass[1] x[1] y[1] z[1] ...
...
id[n-1]  mass[n-1] x[n-1] y[n-1] z[n-1] ...
```
という形式である必要がある。

 -o: 出力ファイルを指定する。出力形式は
 

```
time
n_particles
id[0]  mass[0] x[0] y[0] z[0] nnb[0]  cluster_id[0]
id[1]  mass[1] x[1] y[1] z[1] nnb[1]  cluster_id[1]
...
id[n-1]  mass[n-1] x[n-1] y[n-1] z[n-1] nnb[n-1]  cluster_id[n-1]
```
と、粒子の各行に、その粒子の近傍粒子の数と属するクラスタのid を付加したものになる。


 -r: ネイバー半径 デフォルト 0.01

  -h: ヘルプメッセージを出力して終了する。


# 実行例

```
0
10
0 1  0   0 0
1 1  0.1 0 0
2 1  0.2 0 0
3 1  0.3 0 0
4 1  0.4 0 0
5 1  0.6 0 0
6 1  0.7 0 0
7 1  0.8 0 0
8 1  0.9 0 0
9 1  1.0 0 0
```
がファイル testin にあるとすると
```
fof -i testin -o test.out -r 0.11
```
の実行の結果、 test.out の中身は
```
0.000000e+00
10
0	1	0	0	0	2	0
1	1	0.1	0	0	3	0
2	1	0.2	0	0	3	0
3	1	0.3	0	0	3	0
4	1	0.4	0	0	2	0
5	1	0.6	0	0	2	5
6	1	0.7	0	0	3	5
7	1	0.8	0	0	3	5
8	1	0.9	0	0	3	5
9	1	1	0	0	2	5
```
となるはずである。

# 更新履歴

2022/1/9

* 作ってこれ書いた
