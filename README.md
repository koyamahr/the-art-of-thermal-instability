# The Art of Thermal Instability

熱的不安定性の歴史を、図・コード・計算データからたどる研究アーカイブ。

## 現在地

初期構築版。論文台帳と展示の構成を用意しています。主要6編のPDFをADS・arXivから取得済みです。Field 1965 Figures 1・2・3を式から再計算し、それぞれ15項目・10項目・24項目の検証を通過しました。英語の展示ページ、Pythonコード、CSVと検証記録を用意しています。GitHub Pagesに初版を公開済みです。ライセンス選定は未完了です。

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

取得元・版・チェックサムは `sources/manifest.json` に記録しています。書誌の表題と著者は各PDFの先頭ページで確認済みです。式・図番号とarXiv版／出版版の差異は今後照合します。原論文PDFは自動では公開しません。

ライセンスは未設定です。コード・解説・計算データそれぞれの公開条件は最初の計算公開までに決定します。
