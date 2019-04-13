multiple.py

すべての生物種のゲノムに含まれる「16S r-RNA」という遺伝子の塩基配列を比較し、生物同士の遺伝的な類似度を計算します。
生物種を複数登録することで、その集合全体に対して系統樹を再現します。


"""
Usage:
    1. input: <academic name>
    2. input: <command(arg)>
    
Commands:
    nametable()             Show nametable
    sequence(int)           Show DNA-sequence
    alignment(int, int)     Show alignment result
    cluster()               Show clustering result
    end()                   Stop inputting and start analizing
    save()                  Save to text file(same directory)
    exit()                  Exit this application
"""

input: <academic name>
    ex)input: homo sapiens
    
    生物の学名を入力してください。慣用名(「human」等)で検索できる場合もあります。
    タイプミス・スペルミスは自動で正解を推測してくれます。
    
    入力してEnterを押してもまだ何も出力されません。
    後でまとめて解析するために一時的に配列に登録されます。
    
    一通り登録できたら、「end()」で検索及び解析を開始します。
    
nametable()             Show nametable
    検索済の生物種名が通し番号付きで一覧で表示されます。
    
sequence(int)           Show DNA-sequence
    ex)input: sequence(0)
    
    指定した番号の生物種のDNA(16S r-RNA)配列を表示します。
    
alignment(int, int)     Show alignment result
    ex)input: alignment(0,1)
    
    指定した番号組の生物種間のDNA(16S r-RNA)配列に対して、アラインメント及び類似度を表示します。
    
cluster()               Show clustering result
    解析済の生物種全体に対して、類似度を元にクラスタリングを行いその結果を表示します。
    
save()                  Save to text file(same directory)
    解析済の生物種全体におけるすべての組に関する類似度情報をまとめたものを、テキストファイルに出力します。
    
    