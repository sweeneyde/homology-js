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
    ['globe', [(
        "0: p, q\n" +
        "1: e, f\n" +
        "2: N, S\n"
    ), (
        "e: q - p\n" +
        "f: q - p\n" +
        "N: f - e\n" +
        "S: f - e"
    )]],
    ['sphere_cw', [(
        "0:p; 2:S"
    ), (
        "S:0"
    )]],
    ['torus_simp', [(
        "0: v00,v01,v02, v10,v11,v12, v20,v21,v22\n" +
        "1: a00,a01,a02, a10,a11,a12, a20,a21,a22, b00,b01,b02, b10,b11,b12, b20,b21,b22, c00,c01,c02, c10,c11,c12, c20,c21,c22\n" +
        "2: U00,U01,U02, U10,U11,U12, U20,U21,U22, L00,L01,L02, L10,L11,L12, L20,L21,L22"
    ), (
        "a00:v01-v00; a01:v02-v01; a02:v00-v02\n" +
        "a10:v11-v10; a11:v12-v11; a12:v10-v12\n" +
        "a20:v21-v20; a21:v22-v21; a22:v20-v22\n" +
        "b00:v10-v00; b10:v20-v10; b20:v00-v20\n" +
        "b01:v11-v01; b11:v21-v11; b21:v01-v21\n" +
        "b02:v12-v02; b12:v22-v12; b22:v02-v22\n" +
        "c00:v11-v00; c01:v12-v01; c02:v10-v02\n" +
        "c10:v21-v10; c11:v22-v11; c12:v20-v12\n" +
        "c20:v01-v20; c21:v02-v21; c22:v00-v22\n" +
        "U00:b01-c00+a00; U01:b02-c01+a01; U02:b00-c02+a02\n" +
        "U10:b11-c10+a10; U11:b12-c11+a11; U12:b10-c12+a12\n" +
        "U20:b21-c20+a20; U21:b22-c21+a21; U22:b20-c22+a22\n" +
        "L00:a10-c00+b00; L10:a20-c10+b10; L20:a00-c20+b20\n" +
        "L01:a11-c01+b01; L11:a21-c11+b11; L21:a01-c21+b21\n" +
        "L02:a12-c02+b02; L12:a22-c12+b12; L22:a02-c22+b22"
    )]],
    ['torus_delta', [(
        "0: v\n" +
        "1: a, b, c\n" +
        "2: U, L"
    ), (
        "a: v - v\n" +
        "b: v - v\n" +
        "c: v - v\n" +
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
    ['genus2_delta', [(
        "0: o, x\n" +
        "1: a1,a2,b1,b2, r1,s1,t1,u1, r2,s2,t2,u2\n" +
        "2: R1,S1,T1,U1, R2,S2,T2,U2"
    ), (
        "a1:x-x; b1:x-x; a2:x-x; b2:x-x\n" +
        "r1:x-o; s1:x-o; t1:x-o; u1:x-o\n" +
        "r2:x-o; s2:x-o; t2:x-o; u2:x-o\n" +
        "R1:a1-s1+r1; S1:b1-t1+s1\n" +
        "T1:a1-t1+u1; U1:b1-u1+r2\n" +
        "R2:a2-s2+r2; S2:b2-t2+s2\n" +
        "T2:a2-t2+u2; U2:b2-u2+r1"
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
    ['klein_circle', [(
        "0: vv\n" +
        "1: av, bv, ve\n" +
        "2: Kv, ae, be\n" +
        "3: Ke"
    ), (
        "av:vv-vv; bv:vv-vv; ve:vv-vv\n" +
        "Kv:2av; ae:ve-ve-av+av; be:ve-ve-bv+bv\n" +
        "Ke:2ae+Kv-Kv"
    )]],
    ['klein_klein', [(
        "0: vv\n" +
        "1: av, bv, va, vb\n" +
        "2: Kv, vK, aa, ab, ba, bb\n" +
        "3: aK, bK, Ka, Kb\n" +
        "4: KK"
    ), (
        "av:vv-vv; bv:vv-vv; va:vv-vv; vb:vv-vv\n" +
        "Kv:2av; vK:2va\n" +
        "aa:va-va-av+av; ab:vb-vb-av+av\n" +
        "ba:va-va-bv+bv; bb:vb-vb-bv+bv\n" +
        "aK:vK-vK-2aa; bK:vK-vK-2ba\n" +
        "Ka:2aa+Kv-Kv; Kb:2ab+Kv-Kv\n" +
        "KK:2aK+2Ka"
    )]],
    ['rp2_circle', [(
        "0: vv\n" +
        "1: av, ve\n" +
        "2: Pv, ae\n" +
        "3: Pe"
    ), (
        "av: vv-vv\n" +
        "ve: vv-vv\n" +
        "Pv: 2av\n" +
        "ae: ve-ve-av+av\n" +
        "Pe: 2ae+Pv-Pv"
    )]],
    ['sphere_circle', [(
        "0: vv\n" +
        "1: ve\n" +
        "2: Sv\n" +
        "3: Se"
    ), (
        "ve:vv-vv\n" +
        "Sv:0\n" +
        "Se:Sv-Sv"
    )]],
    ['rp3', [(
        "0: v\n" +
        "1: e\n" +
        "2: f\n" +
        "3: X"
    ), (
        "e: v - v\n" +
        "f: e + e\n" +
        "X: f - f"
    )]],
    ['S3', [(
        "0: v\n" +
        "3: S"
    ), (
        "S: 0"
    )]],
    ['z2_nerve', [(
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
    ['z3_nerve', [(
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

