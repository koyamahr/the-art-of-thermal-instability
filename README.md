# The Art of Thermal Instability

熱的不安定性の歴史を、図・コード・計算データからたどる研究アーカイブ。

## 現在地

初期公開版。論文台帳と展示の構成を用意しています。主要6編のPDFをADS・arXivから取得済みです。Field 1965 Figures 1・2・3を式から再計算し、それぞれ15項目・10項目・24項目の検証を通過しました。FGH/GHF 1969の加熱・金属冷却の計算も追加し、簡略化した低温枝の連立平衡計算も公開し、合計216項目を検証します。Wolfire 1995/2003、GOW 2017と公開コード、化学実装2023、Inoue・Iwasakiの関連10編を英語のResearch pathsに登録しました。反応係数の歴史を独立系譜とし、Takayanagi–Nishimura 1960、Hayakawa–Nishimura–Takayanagi 1961、Hattori–Nakano–Hayashi 1969も追加。相転移フロントの系譜5編（Zel’dovich–Pikel’ner、Penston–Brown、Nagashima 2005/2006、Koyama–Inutsuka 2006 preprint）も登録し全28編。Takayanagi–Nishimura 1960の励起係数34点と有限準位モデルを公開し、原表の冷却率60点との不一致も記録しました。Iwasaki & Inutsuka 2012のFigures 5/7を再現対象に確定し、9例の分類が一致する試作を追加しました。出発分岐・温度範囲の制約とjcritの未確認点を明記し、完全再現とは区別します。FGHの冷温両枝をつなぐS字圧力曲線の試作も追加し、高速電子と二次電離を含む48項目の数値検証を通過しました。水素励起断面積などの仮定を残すため、原図の定量再現は未完了です。英語の展示ページ、Pythonコード、CSVと検証記録を用意しています。GitHub Pagesに初版を公開済みです。ライセンス選定は未完了です。近年の多次元計算は、文献・公開コード／データの案内、取得可能な出力の再解析、小規模検証を基本スコープとし、本計算の再実行は必須にしません。

## ローカルで見る

Python 3で `python3 -m http.server 8000 --directory docs` を実行し、http://localhost:8000 を開きます。サイトにビルドや外部ライブラリは不要です。

## GitHub Pages

公開練習と自動更新の手順は [PUBLISHING.md](PUBLISHING.md) にあります。GitHub Actionsで全図を再計算・検証し、`docs/`だけを公開する構成です。手動の練習モードに加え、mainへのpush後に検証・公開する自動更新を有効にしています。リポジトリは https://github.com/koyamahr/the-art-of-thermal-instability に作成済みです。サイトは https://koyamahr.github.io/the-art-of-thermal-instability/ で公開しています。

## コンテンツ

- `docs/`：公開サイト。`docs/catalog.json` が論文台帳の正本です。
- `ROADMAP.md`：段階的な公開計画と完了条件。
- `CONTRIBUTING.md`：図の再現と計算の記録方針。
- `sources/private/`：作業用PDF。Git管理の対象外です。
- `experiments/`：論文ごとの計算・図生成コードを追加する場所。

取得元・版・チェックサムは `sources/manifest.json` に記録しています。主要6編の表題と著者は各PDFの先頭ページで確認済みです。追加資料の確認範囲は台帳に記録します。式・図番号とarXiv版／出版版の差異は今後照合します。原論文PDFは自動では公開しません。

ライセンスは未設定です。コード・解説・計算データそれぞれの公開条件は最初の計算公開までに決定します。
