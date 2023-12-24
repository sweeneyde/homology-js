let EXAMPLES = new Map([
    ['hollow tetrahedron', [(
        "0: a, b, c, d\n" +
        "1: ab, ac, ad, bc, bd, cd\n" +
        "2: abc, abd, acd, bcd"
    ), (
        "ab: b-a; ac:c-a; ad:d-a; bc:c-b; bd:d-b; cd:d-c\n" +
        "abc:bc-ac+ab; abd:bd-ad+ab; acd:cd-ad+ac; bcd:cd-bd+bc"
    )]],
    ['hollow cube', [(
        "0:AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB\n" +
        "1:AAx,ABx,BAx,BBx,AxA,AxB,BxA,BxB,xAA,xAB,xBA,xBB\n" +
        "2:Axx,Bxx,xAx,xBx,xxA,xxB"
    ), (
        "AAx:AAB-AAA; ABx:ABB-ABA; BAx:BAB-BAA; BBx:BBB-BBA\n" +
        "AxA:ABA-AAA; AxB:ABB-AAB; BxA:BBA-BAA; BxB:BBB-BAB\n" +
        "xAA:BAA-AAA; xAB:BAB-AAB; xBA:BBA-ABA; xBB:BBB-ABB\n" +
        "Axx:AxA+ABx-AxB-AAx; xAx:xAA+BAx-xAB-AAx; xxA:xAA+BxA-xBA-AxA\n" +
        "Bxx:BxA+BBx-BxB-BAx; xBx:xBA+BBx-xBB-ABx; xxB:xAB+BxB-xBB-AxB"
    )]],
    ['hollow octahedron', [(
        "0:Axx,Bxx,xAx,xBx,xxA,xxB\n" +
        "1:AAx,ABx,BAx,BBx,AxA,AxB,BxA,BxB,xAA,xAB,xBA,xBB\n" +
        "2:AAA,AAB,ABA,ABB,BAA,BAB,BBA,BBB"
    ), (
        "AAx:xAx-Axx; ABx:xBx-Axx; BAx:xAx-Bxx; BBx:xBx-Bxx\n" +
        "AxA:xxA-Axx; AxB:xxB-Axx; BxA:xxA-Bxx; BxB:xxB-Bxx\n" +
        "xAA:xxA-xAx; xAB:xxB-xAx; xBA:xxA-xBx; xBB:xxB-xBx\n" +
        "AAA:xAA-AxA+AAx; AAB:xAB-AxB+AAx; ABA:xBA-AxA+ABx\n" +
        "ABB:xBB-AxB+ABx; BAA:xAA-BxA+BAx; BAB:xAB-BxB+BAx\n" +
        "BBA:xBA-BxA+BBx; BBB:xBB-BxB+BBx"
    )]],
    ['sphere_cw', [(
        "0:p; 2:S"
    ), (
        "S:0"
    )]],
    ['torus_simp', [(
        "0: a, b, c, d, e, f, g\n" +
        "1: ab,ac,ad,ae,af,ag,bc,bd,be,bf,bg,cd,ce,cf,cg,de,df,dg,ef,eg,fg\n" +
        "2: abc,abf,acg,ade,adf,aeg,bce,bde,bdg,bfg,cdf,cdg,cef,efg"
    ), (
        "ab:b-a; ac:c-a; ad:d-a; ae:e-a; af:f-a; ag:g-a; bc:c-b\n" +
        "bd:d-b; be:e-b; bf:f-b; bg:g-b; cd:d-c; ce:e-c; cf:f-c\n" +
        "cg:g-c; de:e-d; df:f-d; dg:g-d; ef:f-e; eg:g-e; fg:g-f\n" +
        "abc:bc-ac+ab; abf:bf-af+ab; acg:cg-ag+ac; ade:de-ae+ad\n" +
        "adf:df-af+ad; aeg:eg-ag+ae; bce:ce-be+bc; bde:de-be+bd\n" +
        "bdg:dg-bg+bd; bfg:fg-bg+bf; cdf:df-cf+cd; cdg:dg-cg+cd\n" +
        "cef:ef-cf+ce; efg:fg-eg+ef"
    )]],
    ['torus_delta', [(
        "0: v\n" +
        "1: a, b, c\n" +
        "2: U, L"
    ), (
        "a: 0\n" +
        "b: 0\n" +
        "c: 0\n" +
        "U: b - c + a\n" +
        "L: a - c + b"
    )]],
    ['torus_cw', [(
        "0: v\n" +
        "1: a, b\n" +
        "2: T"
    ), (
        "a:0; b:0;\n" +
        "T: a + b - a - b"
    )]],
    ['rp2_simp', [(
        "0: a,b,c,d,e,f\n" +
        "1: ab,ac,ad,ae,af,bc,bd,be,bf,cd,ce,cf,de,df,ef\n" +
        "2: abc,abf,acd,ade,aef,bce,bde,bdf,cdf,cef"
    ), (
        "ab:b-a; ac:c-a; ad:d-a; ae:e-a; af:f-a; bc:c-b\n" +
        "bd:d-b; be:e-b; bf:f-b; cd:d-c; ce:e-c; cf:f-c\n" +
        "de:e-d; df:f-d; ef:f-e\n" +
        "abc:bc-ac+ab; abf:bf-af+ab; acd:cd-ad+ac; ade:de-ae+ad\n" +
        "aef:ef-af+ae; bce:ce-be+bc; bde:de-be+bd; bdf:df-bf+bd\n" +
        "cdf:df-cf+cd; cef:ef-cf+ce"
    )]],
    ['rp2_delta', [(
        "0: v, w\n" +
        "1: a, b, c\n" +
        "2: U, L"
    ), (
        "a: w-v; b: w-v; c: v-v\n" +
        "U: b - a + c; L: a - b + c"
    )]],
    ['rp2_cw', [(
        "0: v\n" +
        "1: a\n" +
        "2: P"
    ), (
        "P: 2a"
    )]],
    ['T3', [(
        "0: v\n" +
        "1: x, y, z\n" +
        "2: xy, xz, yz\n" +
        "3: xyz"
    ), (
        "x:0; y:0; z:0; xy:0; xz:0; yz:0; xyz:0"
    )]],
    ['genus2', [(
        "0: p\n" +
        "1: a1, b1, a2, b2\n" +
        "2: M"
    ), (
        "M: a1 + b1 - a1 - b1 + a2 + b2 - a2 - b2\n"
    )]],
    ['genus3', [(
        "0: p\n" +
        "1: a1, b1, a2, b2, a3, b3\n" +
        "2: M"
    ), (
        "M: a1 + b1 - a1 - b1 + a2 + b2 - a2 - b2 + a3 + b3 - a3 - b3\n"
    )]],
    ['klein', [(
        "0: v\n" +
        "1: a, b\n" +
        "2: K"
    ), (
        "a: 0\n" +
        "b: 0\n" +
        "K: 2a"
    )]],
    ['z2_bar', [(
        "0:p\n" +
        "1:a\n" +
        "2:aa\n" +
        "3:aaa\n" +
        "4:aaaa\n" +
        "5:aaaaa\n" +
        "6:aaaaaa\n" +
        "7:aaaaaaa\n" +
        "8:aaaaaaaa\n" +
        "9:aaaaaaaaa\n" +
        "10:aaaaaaaaaa\n"
    ), (
        "a:0\n" +
        "aa:a+a\n" +
        "aaa:aa-aa\n" +
        "aaaa:aaa+aaa\n" +
        "aaaaa:aaaa-aaaa\n" +
        "aaaaaa:aaaaa+aaaaa\n" +
        "aaaaaaa:aaaaaa-aaaaaa\n" +
        "aaaaaaaa:aaaaaaa+aaaaaaa\n" +
        "aaaaaaaaa:aaaaaaaa-aaaaaaaa\n" +
        "aaaaaaaaaa:aaaaaaaaa+aaaaaaaaa"
    )]],
    ['z3_bar', [(
        "0: p\n" +
        "1: a, b\n" +
        "2: aa, ab, ba, bb\n" +
        "3: aaa, aab, aba, abb, baa, bab, bba, bbb\n" +
        "4: aaaa, aaab, aaba, aabb, abaa, abab, abba, abbb, baaa, baab, baba, babb, bbaa, bbab, bbba, bbbb"
    ), (
        "aa:a-b+a; ab:b+a; ba:a+b; bb:b-a+b\n" +
        "aaa:aa-ba+ab-aa; aab:ab-bb-aa; aba:ba-ab; abb:bb+aa-ab\n" +
        "baa:aa+bb-ba; bab:ab-ba; bba:ba-aa-bb; bbb:bb-ab+ba-bb\n" +
        "aaaa:aaa-baa+aba-aab+aaa; aaab:aab-bab+abb+aaa; aaba:aba-bba+aab\n" +
        "aabb:abb-bbb-aaa+aab; abaa:baa-abb+aba; abab:bab+aba\n" +
        "abba:bba+aaa+abb; abbb:bbb+aab-aba+abb; baaa:aaa+bba-bab+baa\n" +
        "baab:aab+bbb+baa; baba:aba+bab; babb:abb-baa+bab\n" +
        "bbaa:baa-aaa-bbb+bba; bbab:bab-aab+bba; bbba:bba-aba+baa+bbb\n" +
        "bbbb:bbb-abb+bab-bba+bbb"
    )]],
    ['torus_reduced', [(
        "-1: _\n" +
        "0: v\n" +
        "1: a, b\n" +
        "2: T"
    ), (
        "v: _\n" +
        "a:0; b:0;\n" +
        "T: a + b - a - b"
    )]],
    ['example_1', [(
        "0: p, q\n" +
        "1: x, y, z\n" +
        "2: F"
    ), (
        "x:p-q; y:p-q; z:p-q\n" +
        "F: 2x-2y\n"
    )]],
    ['tietze_result', [(
        "0: v\n" +
        "1: a, b\n" +
        "2: Relation1"
    ), (
        "Relation1: 36b - 28a"
    )]],
]);

