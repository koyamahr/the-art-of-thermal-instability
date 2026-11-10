# The Art of Thermal Instability

KI2002の訂正後の加熱冷却関数から平衡曲線・等圧安定区間を追加。ZP飽和圧力2822.884297 K cm^-3と冷相・不安定枝・温相の温度密度をCSV/JSONに保存。Iwasaki2012実装と照合、27項目を検証。出版版PDFへのアクセスは未成功、独立の正式erratumは未発見と明記。FGHの通常組成は既定範囲の有限探索で飽和圧力を挟めずnull、減耗a/bの暖相拡張では条件付き値を記録（歴史的再現値ではない）。モデルごとの規約・飽和点・未確定理由をdocs/equilibrium-landmarks.jsonに集約。Maxwell1875原典とvan der Waals1873への番外編も追加。

熱的不安定性の歴史を、図・コード・計算データからたどる研究アーカイブ。

## 現在地

初期公開版。論文台帳と展示の構成を用意しています。主要6編のPDFをADS・arXivから取得済みです。Field 1965 Figures 1・2・3を式から再計算し、それぞれ15項目・10項目・24項目の検証を通過しました。FGH/GHF 1969の加熱・金属冷却の計算も追加し、簡略化した低温枝の連立平衡計算も公開し、合計405項目を検証します。Wolfire 1995/2003、GOW 2017と公開コード、化学実装2023、Inoue・Iwasakiの関連10編を英語のResearch pathsに登録しました。反応係数の歴史を独立系譜とし、Takayanagi–Nishimura 1960、Hayakawa–Nishimura–Takayanagi 1961、Hattori–Nakano–Hayashi 1969も追加。相転移フロントの系譜9編（Zel’dovich–Pikel’ner、Penston–Brown、Elphick–Regev–Shaviv 1992、Ferrara–Shchekinov 1993、Shchekinov–Ibáñez 2001、Nagashima 2005/2006、Koyama–Inutsuka 2006 preprint、Kim–Kim 2013）も登録し全47編（教科書1冊を含む）。Hennebelle–Pérault 1999/2000、Burkert–Lin 2000、Audit–Hennebelle 2005、Heitsch et al. 2006、Hennebelle–Audit 2007、Heitsch–Hartmann–Burkert 2008の7編を凝縮・乱流・破砕の系譜として追加。Takayanagi–Nishimura 1960の励起係数34点と有限準位モデルを公開し、原表の冷却率60点との不一致も記録しました。Iwasaki & Inutsuka 2012のFigures 5/7を再現対象に確定し、9例の分類が一致する試作を追加しました。出発分岐・温度範囲の制約とjcritの未確認点を明記し、完全再現とは区別します。FGHの冷温両枝をつなぐS字圧力曲線の試作も追加し、高速電子と二次電離を含む48項目の数値検証を通過しました。水素励起断面積などの仮定を残すため、原図の定量再現は未完了です。英語の展示ページ、Pythonコード、CSVと検証記録を用意しています。GitHub Pagesに初版を公開済みです。ライセンス選定は未完了です。近年の多次元計算は、文献・公開コード／データの案内、取得可能な出力の再解析、小規模検証を基本スコープとし、本計算の再実行は必須にしません。

## ローカルで見る

Python 3の標準ライブラリだけで公開用HTMLを生成できます。

```sh
python3 scripts/build_site.py
python3 scripts/check_site.py --directory dist
python3 -m http.server 8000 --bind 127.0.0.1 --directory dist
```

http://localhost:8000 を開きます。既存18ページに加え、文献47件の個別紹介を `papers/<id>/` に生成します。一覧はJavaScriptなしでも読め、スクリプトはトピック絞り込みを補助します。既存の図・計算ページと `#field-1965` などのアンカーは引き継ぎます。全65ページに固有のタイトル・説明文、canonical、Open Graph情報を付け、`sitemap.xml` を生成します。英語の既存資料に対応する構成で、日本語訳は追加していません。

## GitHub Pages

公開練習と自動更新の手順は [PUBLISHING.md](PUBLISHING.md) にあります。GitHub Actionsで全図を再計算・検証し、公開用 `docs/` から静的HTMLを生成し、検証済みの `dist/` だけを公開する構成です。手動の練習モードに加え、mainへのpush後に検証・公開する自動更新を有効にしています。リポジトリは https://github.com/koyamahr/the-art-of-thermal-instability に作成済みです。サイトは https://koyamahr.github.io/the-art-of-thermal-instability/ で公開しています。

## コンテンツ

- `docs/`：公開サイト。`docs/catalog.json` が論文台帳の正本です。
- `ROADMAP.md`：段階的な公開計画と完了条件。
- `CONTRIBUTING.md`：図の再現と計算の記録方針。
- `sources/private/`：作業用PDF。Git管理の対象外です。
- `experiments/`：論文ごとの計算・図生成コードを追加する場所。

取得元・版・チェックサムは `sources/manifest.json` に記録しています。主要6編の表題と著者は各PDFの先頭ページで確認済みです。追加資料の確認範囲は台帳に記録します。式・図番号とarXiv版／出版版の差異は今後照合します。原論文PDFは自動では公開しません。

自作の実装コードはMIT。解説・図・計算データ・動画はMITの対象外で、個別条件または権利留保です。詳細はLICENSING.mdを参照してください。

追加基礎文献：Parker 1953、Spitzer–Tomasko 1968、Dalgarno–McCray 1972、Balbus 1986、Begelman–McKee 1990、Bakes–Tielens 1994。新規の非線形TIシミュレーション展示では熱伝導と対応する解像度の扱いを選定時に確認する。解析理論・加熱冷却係数の原典は別の役割として整理する。

Kim & Kim (2013, ApJ779,48)はユーザーの追加指定により登録。熱伝導を含む蒸発フロントの線形安定性・非線形飽和の研究として、Inoue2006とIwasaki2012への参照を追加。

FGHは原図96点との比較・係数監査・80項目の数値チェックを追加。通常組成の限界圧力は1516.7対1800Kcm^-3で約15.7%低く、遷移位置のずれにより局所差はさらに大きい。原因を特定できないため、ユーザーの意向に従い、厳密な歴史的係数の回収は新しい根拠が得られるまで保留。比較図・生ピクセル座標・較正・CSV・加熱冷却内訳・感度診断・再開条件を保存する。

Schwarz, McCray & Stein (1972)のFigure2を最初の再現対象とし、冷却・再結合背景と熱伝導を含む線形試作を追加。Figure1の39読み取り点を冷却入力とし、Figure2の24点と比較。入力下限13.929Kに達する1.679Myrで停止し、後半と非線形Figures3–6は未再現と明記。独立積分法・元方程式のヤコビアン等31項目を検証。

Schwarz1972の非線形ラグランジュ試作も公開。リーマンソルバーを用いず、互換圧力仕事・圧縮時人工粘性・熱伝導・再結合を実装。64/96/128メッシュと補助512セルを比較し、表示は96セル。33項目の数値自己検証を通過する一方、細い伝導層の十分な解像は未確認と明示。0.985Myrで密度差約187、最低温度12.7K、密度FWHM0.0078pc。原論文との時刻・温度・電離度の違い、冷却入力外挿を記録。


## Licensing / ライセンス

Original implementation code: **MIT** ([LICENSE](LICENSE)). Text, figures, photographs, book material and data are **excluded** and follow separate terms; see [LICENSING.md](LICENSING.md) and [third-party notices](THIRD_PARTY_NOTICES.md). Public guide: https://koyamahr.github.io/rights/ .
