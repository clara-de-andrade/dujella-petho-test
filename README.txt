python/
    main.py:
        Programa para gerar cache.pkl e reduzir (usando a library mpmath) as
        cotas superiores para os inteiros positivos m, n, s, k satisfazendo
            
            F(k, m) ** s + F(k, m + 1) ** s = F(k, n)
            
        quando 3 <= m <= 1394 e 3 <= k <= min(m, log s), a partir do Teorema de
        Dujella-Petho. (Ver comentários no código fonte.)

    test.py:
        Funções para calcular os valores necessários em python/main.py. (Ver
        comentários no código fonte.)

    base.py:
        Funções miscelâneas usadas em test.py.

    cache.py:
        Funções para armazenar os valores calculados em cache.pkl (usando a
        library pickle).

    cache.pkl:
        Arquivo binário com as cotas calculadas e vários valores importantes
        associados. Execute

        ```
        with open("cache.pkl", "rb") as f:
            cache = pickle.load(f)
        ```

        em uma shell do Python para acessar o conteúdo do arquivo como um
        dicionário `cache`.

    demo.mp4:
        Demonstração do programa main.py em operação. (47min12s)