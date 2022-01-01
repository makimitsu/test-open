# test-open

### testの変更
### branch2を作成


# 簡単なやりかた
## アカウントをつくる
https://qiita.com/ayatokura/items/9eabb7ae20752e6dc79d

## Github Desktopのインストール
https://docs.github.com/ja/desktop/installing-and-configuring-github-desktop/installing-and-authenticating-to-github-desktop/installing-github-desktop
## 編集したいレポジトリを用意する
### repositoryを作成（新規作成、または、ローカルからオンライン）
*DesktopアプリからはFile> New repository（新規作成)　か　Add local repository(ローカルにあるフォルダをレポジトリにする)

### repositoryをクローン（オンラインからローカルで編集できるようにする）
*DesktopアプリからはFile> Clone repository > URL　でオンラインのレポジトリをクローンしてくることができる
*編集するときは権限が必要なので権限のある人からinvitationを送ってもらう

## repositoryを編集する
参考⇒　https://docs.github.com/ja/desktop/contributing-and-collaborating-using-github-desktop

各レポジトリはmainが本体で、各自がそれをクローン（フォルダのコピー）して編集する。この時本体のmainをみんなで編集してごちゃごちゃにならないように、mainからブランチ（編集する際の枝分かれ）を作ってそこで編集をする。編集したブランチを最終的にmainにマージ（合体）させることで本体のmainのコードを発展させていくイメージ。ブランチは始めはmainから枝分かれするが、branchの中でも編集を分けたいときはブランチからさらに枝分かれしたブランチを切ることができる。（２.参照）
1. ブランチを切る。create new branch　⇒　好きな名前のブランチを作成（例：main>branch1）。これでmainもしくはひとつ前のブランチを直接編集せずに枝分かれの先のみ編集できるようになる。
2. ブランチを編集する。
   1. ブランチ内でさらに枝分かれしたいときはブランチからさらにブランチを切れる。（例：mainから枝分かれしたbranch1を見て、このブランチもっとこうしたほうがいいんじゃないみたいに思ったら、branch1からbranch1_1みたいに新たにブランチを切る。branch1_1で編集した内容で大丈夫になったらbranch1にマージ。branch1の編集がこれでOKになったらbranch1をmainにマージという感じ）
3. Push Origin。ローカルで編集したファイルをオンラインのGithubのレポジトリに反映させる。
4. Pull request。反映された更新をmain（もしくはひとつ前のbranch）にマージできるかどうかリクエストする。
5. merge。pull requestでコンフリクトが起きなければ編集したコードをmain（もしくはひとつ前のbranch）にマージさせる。マージがすんだらbranchとマージ先のmain（もしくはひとつ前のbranch）が同じものになり、branchが削除可能になる（そのブランチでの編集の終わり）。
## 誰かが編集したrepositoryをbranchにマージする
誰かがmainにマージすると、その前のバージョンから切ったブランチにはあたらしい更新が反映されていない。編集中のブランチにあたらしい更新をマージする方法。
1. GitHubのwebページに行く
2. Pull request。　branchからmainのときと反対に　main（もしくはひとつ前のブランチ）>branchのpull requestをする
3. merge。コンフリクトがなければマージ。
4. Pull Origin。オンラインのGithubのレポジトリをローカルのファイルに反映させる。これでローカルのファイルが更新されたmain（もしくはひとつ前のブランチ）のものと同じになった

