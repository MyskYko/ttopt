/* Generated by Yosys 0.9 (git sha1 1979e0b) */

module dct(\x[0] , \x[1] , \x[2] , \y[0] , \y[1] , \y[2] , \u[0] , \u[1] , \u[2] , \v[0] , \v[1] , \v[2] , \coef[0] , \coef[1] , \coef[2] , \coef[3] , \coef[4] , \coef[5] , \coef[6] , \coef[7] , \coef[8] , \coef[9] , \coef[10] , \coef[11] , \coef[12] , \coef[13] , \coef[14] , \coef[15] , \coef[16] , \coef[17] , \coef[18] , \coef[19] , \coef[20] , \coef[21] , \coef[22] , \coef[23] , \coef[24] , \coef[25] , \coef[26] , \coef[27] , \coef[28] , \coef[29] , \coef[30] , \coef[31] );
  output \coef[0] ;
  output \coef[10] ;
  output \coef[11] ;
  output \coef[12] ;
  output \coef[13] ;
  output \coef[14] ;
  output \coef[15] ;
  output \coef[16] ;
  output \coef[17] ;
  output \coef[18] ;
  output \coef[19] ;
  output \coef[1] ;
  output \coef[20] ;
  output \coef[21] ;
  output \coef[22] ;
  output \coef[23] ;
  output \coef[24] ;
  output \coef[25] ;
  output \coef[26] ;
  output \coef[27] ;
  output \coef[28] ;
  output \coef[29] ;
  output \coef[2] ;
  output \coef[30] ;
  output \coef[31] ;
  output \coef[3] ;
  output \coef[4] ;
  output \coef[5] ;
  output \coef[6] ;
  output \coef[7] ;
  output \coef[8] ;
  output \coef[9] ;
  input \u[0] ;
  input \u[1] ;
  input \u[2] ;
  input \v[0] ;
  input \v[1] ;
  input \v[2] ;
  input \x[0] ;
  input \x[1] ;
  input \x[2] ;
  input \y[0] ;
  input \y[1] ;
  input \y[2] ;
  assign \coef[0]  = 4096'h3c3cffffffff3c3c0000ff6666ff00005a5affffffff5a5aff00ff0000ff00ffa5a5ffffffffa5a50000ff9999ff0000c3c3ffffffffc3c3ff00ff0000ff00ff243c3c24243c3c249900009999000099185a5a18185a5a18ff0000ffff0000ff81a5a58181a5a581660000666600006642c3c34242c3c342ff0000ffff0000ff3cff3cffff3cff3c006600ffff0066005aff5affff5aff5a0000ffffffff0000a5ffa5ffffa5ffa5009900ffff009900c3ffc3ffffc3ffc30000ffffffff0000a5a5a5a5a5a5a5a599999999999999993c3c3c3c3c3c3c3c007f405f50575455c3c3c3c3c3c3c3c366666666666666665a5a5a5a5a5a5a5a0000000000000000ff3cff3c3cff3cffff006600006600ffff5aff5a5aff5affffff00000000ffffffa5ffa5a5ffa5ffff009900009900ffffc3ffc3c3ffc3ffffff00000000ffff3c24243c3c24243c00999900009999005a18185a5a18185a00ffff0000ffff00a58181a5a58181a50066660000666600c34242c3c34242c300ffff0000ffff00ffff3c3c3c3cffff66ff00000000ff66ffff5a5a5a5affff00ff00ffff00ff00ffffa5a5a5a5ffff99ff00000000ff99ffffc3c3c3c3ffff00ff00ffff00ff00a5a5a5a5a5a5a5a599999999999999993c3c3c3c3c3c3c3c0000000000000000c3c3c3c3c3c3c3c366666666666666665a5a5a5a5a5a5a5a0000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[9]  = 4096'h55d40fd42bf02baaa55ac35aa53ca55a4d698e69967196b2666699669966999971334d33ccb2cc8ec33c5a3cc3a5c33c0f17aa17e855e8f00000ff00ff00ffffaa718e55558e71aa3c5aa5c3c3a55a3cb255aa4d4daa55b299996666666699998ef00f71710ff08ea53cc35a5ac33ca5f04db20f0fb24df0ffff00000000ffff2bd4550ff0aa2bd4a55aa5c33c5aa55a96694d8e71b296699966669966999966cc33714db28ecc33c33cc35aa53cc33ce8170faa55f0e817ff0000ff00ffff002bd4d42b2bd4d42bc33c3cc3c33c3cc396696996966969960066401910464411cc3333cccc3333cc5aa5a55a5aa5a55ae81717e8e81717e800000000000000000faad4d42b2b55f0c35a5a5aa5a5a53c8eb2696996964d7199996666999966664d8e3333cccc71b25a3c3c3cc3c3c3a5aaf01717e8e80f55ffff0000ffff00008eaa55717155aa8ea53cc35a5ac33ca5aab24d55554db2aa66996699996699660f8e71f0f0718e0fc3a55a3c3c5aa5c3b2f00f4d4d0ff0b200ff00ffff00ff002b0f2b55aad4f0d4a5c3a5a55a5a3c5a968e964db26971699999996699666666cc4dcc718e33b233c35ac3c33c3ca53ce8aae80ff0175517ffffff00ff0000002b2b2b2b2b2b2b2bc3c3c3c3c3c3c3c396969696969696960000000000000000cccccccccccccccc5a5a5a5a5a5a5a5ae8e8e8e8e8e8e8e80000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[10]  = 4096'hb271d4e8172b8e4dc3a5a53cc35a5a3cf0556933cc96aa0f9966666699999966aaf0339669cc0f555ac3c3a55a3c3ca5714d17d42be8b28eff000000ffffff00e88e711717718ee83c5aa5c3c3a55a3c33aa55cccc55aa336666999999996666960ff06969f00f96a53cc35a5ac33ca5d4b24d2b2b4db2d40000ffffffff00008ee8b2d42b4d17715a3cc3a55a3cc3a5aa33f069960fcc5599669966996699660f96aa33cc5569f03ca55ac33ca55ac3b2d47117e88e2b4dff00ff00ff00ff008e71718e8e71718e3cc3c33c3cc3c33caa5555aaaa5555aa00664019104644110ff0f00f0ff0f00fa55a5aa5a55a5aa5b24d4db2b24d4db20000000000000000d44de8718e17b22ba53c3ca55ac3c35a690f3355aaccf0966666666699999999335596f00f69aaccc3a5a5c33c5a5a3c178ed44db22b71e800000000ffffffff71e8178e8e17e871a53cc35a5ac33ca55533ccaaaacc33559966996666996699f096690f0f6996f0c3a55a3c3c5aa5c34dd42bb2b22bd44dff00ff0000ff00ff17d48eb24d712be8c3a55ac33ca55a3ccc69aaf00f559633996699996666996669330faa55f0cc965ac33c5aa5c33ca52b17b2718e4de8d4ff00ffff0000ff008e8e8e8e8e8e8e8e3c3c3c3c3c3c3c3caaaaaaaaaaaaaaaa00000000000000000f0f0f0f0f0f0f0fa5a5a5a5a5a5a5a5b2b2b2b2b2b2b2b20000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[11]  = 4096'h33558ed42b71aaccc3c33ca55ac33c3cd44daa699655b22b6699996699666699e8710f33ccf08e175a5aa5c33c5aa5a5690fb217e84df09600ffff00ff0000ffd4cc332b2b33ccd43c5aa5c3c3a55a3c692bd49696d42b6999669966669966993317e8cccce81733a53cc35a5ac33ca5179669e8e8699617ff00ff0000ff00ffaad4338e71cc2b553ca5c33cc33c5ac3b269d4aa552b964d66666699669999998e33e80ff017cc71a5c35aa55aa53c5af01769b24d96e80f000000ff00ffffff6996966969969669a55a5aa5a55a5aa517e8e81717e8e8170066401910464411d42b2bd4d42b2bd4c33c3cc3c33c3cc3cc3333cccc3333cc00000000000000008eccd455aa2b33713c3ca5c33c5ac3c3aa2b694db296d45599996699669966660f1733718ecce8f0a5a5c35aa53c5a5ab296170ff0e8694dffff00ff00ff000033d42bcccc2bd433a53cc35a5ac33ca5d469962b2b9669d49999666666669999e833cc1717cc33e8c3a55a3c3c5aa5c36917e89696e81769ffff00000000ffff2b8eaa33cc5571d45a3c3cc33cc3c3a596aab2d42b4d55699999666699996666cc0f8ee81771f0333ca5a55aa55a5ac3e8b2f069960f4d17ffff0000ffff00006969696969696969a5a5a5a5a5a5a5a517171717171717170000000000000000d4d4d4d4d4d4d4d4c3c3c3c3c3c3c3c3cccccccccccccccc0000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[12]  = 4096'hccaa712bd48e55333ca55a3cc3a55ac32bb2559669aa4dd49966996699669966178ef0cc330f71e8a5c33ca55ac33c5a96f04de817b20f69ff00ff00ff00ff004d2bd4b2b2d42b4d3c5aa5c3c3a55a3c0f9669f0f069960f996699666699669955cc33aaaa33cc55a53cc35a5ac33ca58ee817717117e88eff00ff0000ff00ff552bcc718e33d4aa5a3c3c5aa5c3c3a54d962b55aad469b2996699996666996671cc17f00fe8338e3ca5a53cc35a5ac30fe8964db26917f0ff00ffff0000ff00aa5555aaaa5555aaa55a5aa5a55a5aa5b24d4db2b24d4db200664019104644118e71718e8e71718ec33c3cc3c33c3cc3f00f0ff0f00f0ff0000000000000000071332baa55d4cc8e5ac33ca55ac33ca555d496b24d692baa9966666699999966f0e8cc8e7133170f3c5aa5c33c5aa5c34d69e8f00f1796b2ff000000ffffff00d44db22b2bb24dd4a53cc35a5ac33ca5690ff09696f00f6999996666666699993355aaccccaa5533c3a55a3c3c5aa5c3178e71e8e8718e17ffff00000000ffffd47155cc33aa8e2bc35a5a3cc3a5a53c69554d2bd4b2aa96999999996666666633f07117e88e0fcc5a3c3ca55ac3c3a5174d0f9669f0b2e8ffffffff00000000aaaaaaaaaaaaaaaaa5a5a5a5a5a5a5a5b2b2b2b2b2b2b2b200000000000000008e8e8e8e8e8e8e8ec3c3c3c3c3c3c3c3f0f0f0f0f0f0f0f00000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[13]  = 4096'h698eb233cc4d71963c5a3ca55ac3a5c317aaf0d42b0f55e86699999966666699d40faae81755f02ba53ca5c33c5ac35accb27169968e4d3300ffffff000000ff174db2e8e8b24d17c3a55a3c3c5aa5c3cc0ff03333f00fcc99996666666699996955aa9696aa55695ac33ca5a53cc35a2b8e71d4d4718e2bffff00000000ffff713369b24d96cc8ea5a53c3cc3c35a5a55d417f00fe82baa6699669966996699f0e8d4aa552b170fc3c3a5a55a5a3c3c4d69cc718e3396b200ff00ff00ff00ff718e8e71718e8e71c33c3cc3c33c3cc355aaaa5555aaaa550066401910464411f00f0ff0f00f0ff05aa5a55a5aa5a55a4db2b24d4db2b24d0000000000000000b296338e71cc694d3cc3a55aa55a3cc3f0e8d4aa552b170f9999999966666666aa2be80ff017d455a55ac33cc33ca55a713369b24d96cc8effffffff00000000b217e84d4de817b25ac33ca5a53cc35af0cc330f0f33ccf06699669999669966aa699655559669aa3c5aa5c3c3a55a3c712bd48e8ed42b7100ff00ffff00ff00ccb27169968e4d335a3ca53cc35ac3a52bf05517e8aa0fd4669966669999669917aaf0d42b0f55e83ca5c3a55a3c5ac396714dcc33b28e6900ff0000ffff00ff7171717171717171c3c3c3c3c3c3c3c355555555555555550000000000000000f0f0f0f0f0f0f0f05a5a5a5a5a5a5a5a4d4d4d4d4d4d4d4d0000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[14]  = 4096'h0f0f3333ccccf0f0a5c3c3a55a3c3c5a8e8ed4d42b2b717199996666999966664d4de8e81717b2b2c35a5ac33ca5a53caaaa696996965555ffff0000ffff0000f069960f0f9669f0c3a55a3c3c5aa5c37117e88e8ee817719966996666996699b2d42b4d4d2bd4b25ac33ca5a53cc35a55cc33aaaa33cc55ff00ff0000ff00fff0330f33ccf0cc0f3ca5a5c33c5a5ac371d48ed42b712b8e6666996699669999b2e84de817b2174da5c3c35aa53c3c5a5569aa69965596aa0000ff00ff00ffffcc3333cccc3333cca55a5aa5a55a5aa52bd4d42b2bd4d42b006640191046441117e8e81717e8e817c33c3cc3c33c3cc39669699696696996000000000000000033f0330ff0cc0fccc35aa5c33c5aa53cd471d48e712b8e2b6666669966999999e8b2e84db2174d175a3cc35aa53cc3a5695569aa5596aa96000000ff00ffffff96f00f69690ff0965ac33ca5a53cc35ae8718e17178e71e899996666666699992bb24dd4d44db22b3c5aa5c3c3a55a3c3355aaccccaa5533ffff00000000ffffcc33f00ff00fcc335ac33ca55ac33ca52bd4718e718e2bd4996666996699996617e8b24db24d17e83c5aa5c33c5aa5c3966955aa55aa9669ff0000ff00ffff00cccccccccccccccca5a5a5a5a5a5a5a52b2b2b2b2b2b2b2b00000000000000001717171717171717c3c3c3c3c3c3c3c396969696969696960000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[15]  = 4096'haa33f069960fcc55a53cc33cc33cc35ab2d47117e88e2b4d99666699669999668ee8b2d42b4d1771c3a55aa55aa55a3cf06955cc33aa960fff0000ff00ffff00aa2bd45555d42baa3c5aa5c3c3a55a3cb296694d4d6996b266996699996699668ecc33717133cc8ea53cc35a5ac33ca5f0e8170f0f17e8f000ff00ffff00ff00cc69aaf00f559633c33ca5c33c5ac33c2b17b2718e4de8d4999999669966666617d48eb24d712be85aa5c35aa53c5aa596ccf055aa0f3369ffffff00ff00000096696996966969965aa5a55a5aa5a55ae81717e8e81717e800664019104644112bd4d42b2bd4d42b3cc3c33c3cc3c33c33cccc3333cccc330000000000000000f0556933cc96aa0fc35a3c3cc3c3a53c714d17d42be8b28e6666996699669999b271d4e8172b8e4d5a3ca5a55a5ac3a5550fcc699633f0aa0000ff00ff00ffffd4aa552b2b55aad4a53cc35a5ac33ca569b24d96964db2696666999999996666338e71cccc718e33c3a55a3c3c5aa5c317f00fe8e80ff0170000ffffffff000096f0ccaa55330f69c3c3c3a55a3c3c3ce8712bb24dd48e1766669999666699992bb2178e71e84dd45a5a5ac33ca5a5a5335596f00f69aacc0000ffff0000ffff96969696969696965a5a5a5a5a5a5a5ae8e8e8e8e8e8e8e800000000000000002b2b2b2b2b2b2b2b3c3c3c3c3c3c3c3c33333333333333330000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[16]  = 4096'hb271cccc33338e4d5aa55aa55aa55aa5f0552b2bd4d4aa0f9966996699669966aaf01717e8e80f553cc33cc33cc33cc3714d96966969b28eff00ff00ff00ff0055aa55aaaa55aa55c3a55a3c3c5aa5c34db24db2b24db24d9966996666996699718e718e8e718e715ac33ca5a53cc35a0ff00ff0f00ff00fff00ff0000ff00ff8eccb2cc334d33715aa55a5aa5a55aa5aa2bf02bd40fd45599669999666699660f17aa17e855e8f03cc33c3cc3c33cc3b2967196698e694dff00ffff0000ff00aa5555aaaa5555aaa55a5aa5a55a5aa5b24d4db2b24d4db200664019104644118e71718e8e71718ec33c3cc3c33c3cc3f00f0ff0f00f0ff00000000000000000cc4dcc718e33b2335aa5a5a55a5a5aa52b0f2b55aad4f0d49966666699999966175517f00fe8aae83cc3c3c33c3c3cc3968e964db2697169ff000000ffffff005555aaaaaaaa55555ac33ca5a53cc35a4d4db2b2b2b24d4d999966666666999971718e8e8e8e71713c5aa5c3c3a55a3c0f0ff0f0f0f00f0fffff00000000ffff33cc8eb24d7133cc5a5a5a5aa5a5a5a5d42baaf00f55d42b9999999966666666e8170faa55f0e8173c3c3c3cc3c3c3c36996b2718e4d6996ffffffff00000000aaaaaaaaaaaaaaaaa5a5a5a5a5a5a5a5b2b2b2b2b2b2b2b200000000000000008e8e8e8e8e8e8e8ec3c3c3c3c3c3c3c3f0f0f0f0f0f0f0f00000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[17]  = 4096'h8e337169968ecc713cc3a55aa55a3cc3aad45517e8aa2b5599669999666699660fe8f0d42b0f17f0a55ac33cc33ca55ab2694dcc33b2964dff00ffff0000ff006955aa9696aa5569c3a55a3c3c5aa5c3174db2e8e8b24d176666999999996666d4718e2b2b8e71d45ac33ca5a53cc35acc0ff03333f00fcc0000ffffffff0000cc698e718e7196333c5a3ca55ac3a5c32b17aa55aa55e8d4999999996666666617d40ff00ff02be8a53ca5c33c5ac35a96ccb24db24d3369ffffffff00000000b24d4db2b24d4db23cc3c33c3cc3c33cf00f0ff0f00f0ff00066401910464411aa5555aaaa5555aaa55a5aa5a55a5aa5718e8e71718e8e71000000000000000071716933cc968e8ea5c35ac33ca53c5a555517d42be8aaaa9966996699669966f0f0d4e8172b0f0fc35a3c5aa5c3a53c4d4dcc699633b2b2ff00ff00ff00ff00aa699655559669aa5ac33ca5a53cc35ab217e84d4de817b299669966669966998ed42b71712bd48e3c5aa5c3c3a55a3cf0cc330f0f33ccf0ff00ff0000ff00ff9671cc8e71338e69a5a53c3cc3c35a5ae8552baa55d4aa1766999999666666992bf0170ff0e80fd4c3c3a5a55a5a3c3c334d96b24d69b2cc00ffffff000000ffb2b2b2b2b2b2b2b23c3c3c3c3c3c3c3cf0f0f0f0f0f0f0f00000000000000000aaaaaaaaaaaaaaaaa5a5a5a5a5a5a5a571717171717171710000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[18]  = 4096'h8e337169968ecc713c5a3ca55ac3a5c3aad45517e8aa2b5566999966996666990fe8f0d42b0f17f0a53ca5c33c5ac35ab2694dcc33b2964d00ffff00ff0000ff174db2e8e8b24d17c3a55a3c3c5aa5c3cc0ff03333f00fcc66669999999966666955aa9696aa55695ac33ca5a53cc35a2b8e71d4d4718e2b0000ffffffff0000cc698e718e719633a5a53c3cc3c35a5a2b17aa55aa55e8d4666666996699999917d40ff00ff02be8c3c3a5a55a5a3c3c96ccb24db24d3369000000ff00ffffff69969669699696693cc3c33c3cc3c33c17e8e81717e8e8170066401910464411d42b2bd4d42b2bd4a55a5aa5a55a5aa5cc3333cccc3333cc000000000000000071716933cc968e8e3cc3a55aa55a3cc3555517d42be8aaaa9999669966996666f0f0d4e8172b0f0fa55ac33cc33ca55a4d4dcc699633b2b2ffff00ff00ff0000b217e84d4de817b25ac33ca5a53cc35af0cc330f0f33ccf09966996666996699aa699655559669aa3c5aa5c3c3a55a3c712bd48e8ed42b71ff00ff0000ff00ff9671cc8e71338e695a3ca53cc35ac3a5e8552baa55d4aa1799996666999966662bf0170ff0e80fd43ca5c3a55a3c5ac3334d96b24d69b2ccffff0000ffff000069696969696969693c3c3c3c3c3c3c3c17171717171717170000000000000000d4d4d4d4d4d4d4d4a5a5a5a5a5a5a5a5cccccccccccccccc0000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[1]  = 4096'h280cb233cc4d3014ff99a542245a99ff120af0d42b0f5048660066ffff9900998405aae81755a021ff66c318813c66ffc0827169968e4103000000ffffff00ffe395a9c7c7a995e37effffe7e7ffff7eb56c36adad366cb599ffff6666ffff99da63c65b5bc663dabdffffdbdbffffbd7c1bd83e3ed81b7cffffff0000ffffff303328b24d14cc0c9942ffa55aff249950d412f00f482b0a00ff66669999ff00a0e884aa552117056618ffc33cff81664169c0718e03968200ff0000ffffff001db8b81d1db8b81de77e7ee7e77e7ee74e72724e4e72724e006640191046441165a6a66565a6a665dbbdbddbdbbdbddb8bd1d18b8bd1d18b0000000000000000b214330c30cc284da5ff42999924ff5af048d40a502b120f6699ff0000ff6699aa21e805a0178455c3ff18666681ff3c710369824196c08e00ffff0000ff00ffa9e3c79595c7e3a9ff7ee7ffffe77eff36b5ad6c6cadb536ff9966ffff6699ffc6da5b63635bdac6ffbddbffffdbbdffd87c3e1b1b3e7cd8ffff00ffff00ffffccb23028140c4d3324a599ffff995a422bf05012480a0fd4ff660066990099ff17aaa084210555e881c366ffff663c18967141c003828e69ff000000ff00ffff1d1d1d1d1d1d1d1de7e7e7e7e7e7e7e74e4e4e4e4e4e4e4e00000000000000006565656565656565dbdbdbdbdbdbdbdb8b8b8b8b8b8b8b8b0000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[19]  = 4096'h8e337169968ecc71c3a55a3cc3a55a3caad45517e8aa2b5599669966996699660fe8f0d42b0f17f05ac33ca55ac33ca5b2694dcc33b2964dff00ff00ff00ff00ccaa55333355aaccc3a55a3c3c5aa5c32bb24dd4d44db22b6666999999996666178e71e8e8718e175ac33ca5a53cc35a96f00f69690ff0960000ffffffff0000cc698e718e7196335a3cc35aa53cc3a52b17aa55aa55e8d4996699996666996617d40ff00ff02be83ca55a3cc3a55ac396ccb24db24d3369ff00ffff0000ff00aa5555aaaa5555aa3cc3c33c3cc3c33cb24d4db2b24d4db200664019104644118e71718e8e71718ea55a5aa5a55a5aa5f00f0ff0f00f0ff0000000000000000071716933cc968e8e5a3c3ca55ac3c3a5555517d42be8aaaa9966666699999966f0f0d4e8172b0f0f3ca5a5c33c5a5ac34d4dcc699633b2b2ff000000ffffff0055cc33aaaa33cc555ac33ca5a53cc35a4d2bd4b2b2d42b4d99669966669966997117e88e8ee817713c5aa5c3c3a55a3c0f9669f0f069960fff00ff0000ff00ff9671cc8e71338e69c35a5ac33ca5a53ce8552baa55d4aa1799999999666666662bf0170ff0e80fd45a3c3c5aa5c3c3a5334d96b24d69b2ccffffffff00000000aaaaaaaaaaaaaaaa3c3c3c3c3c3c3c3cb2b2b2b2b2b2b2b200000000000000008e8e8e8e8e8e8e8ea5a5a5a5a5a5a5a5f0f0f0f0f0f0f0f00000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[20]  = 4096'h4d962b55aad469b2a5a53c3cc3c35a5a0fe8964db26917f06666996699669999552bcc718e33d4aac3c3a5a55a5a3c3c8e33e80ff017cc710000ff00ff00ffffcc0ff03333f00fccc3a55a3c3c5aa5c32b8e71d4d4718e2b9999666666669999174db2e8e8b24d175ac33ca5a53cc35a96aa55696955aa96ffff00000000ffff69554d2bd4b2aa965a3ca53cc35ac3a5174d0f9669f0b2e89966669966999966d47155cc33aa8e2b3ca5c3a55a3c5ac3cc0f8ee81771f033ff0000ff00ffff002bd4d42b2bd4d42bc33c3cc3c33c3cc396696996966969960066401910464411cc3333cccc3333cc5aa5a55a5aa5a55ae81717e8e81717e800000000000000002bb2559669aa4dd43c5a3ca55ac3a5c396f04de817b20f699999666699996666ccaa712bd48e5533a53ca5c33c5ac35ae8710f33ccf08e17ffff0000ffff0000f0cc330f0f33ccf05ac33ca5a53cc35a712bd48e8ed42b716699669999669966b217e84d4de817b23c5aa5c3c3a55a3c559669aaaa69965500ff00ffff00ff00aa2b694db296d455c33c5aa55aa5c33cb296170ff0e8694d99999966996666668eccd455aa2b33715aa53cc33cc35aa5f0e8cc8e7133170fffffff00ff0000002b2b2b2b2b2b2b2bc3c3c3c3c3c3c3c396969696969696960000000000000000cccccccccccccccc5a5a5a5a5a5a5a5ae8e8e8e8e8e8e8e80000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[21]  = 4096'haa2be80ff017d455a55ac33cc33ca55ab296338e71cc694d99996699669966668ecc964db2693371c33c5aa55aa5c33cf0e8d4aa552b170fffff00ff00ff0000aa699655559669aa3c5aa5c3c3a55a3cb217e84d4de817b266996699996699668ed42b71712bd48ea53cc35a5ac33ca5f0cc330f0f33ccf000ff00ffff00ff00d40faae81755f02ba53ca5c33c5ac35a698eb233cc4d71966699996699666699334d8e966971b2ccc3a5c35aa53c5a3c17aaf0d42b0f55e800ffff00ff0000ffd42b2bd4d42b2bd45aa5a55a5aa5a55a6996966969969669006640191046441133cccc3333cccc333cc3c33c3cc3c33c17e8e81717e8e8170000000000000000e8550f2bd4f0aa17c35a3c5aa5c3a53c334d8e966971b2cc666699996666999996714dcc33b28e695a3ca53cc35ac3a5d40faae81755f02b0000ffff0000ffff96aa55696955aa96a53cc35a5ac33ca5e8b24d17174db2e866669999999966662b8e71d4d4718e2bc3a55a3c3c5aa5c333f00fcccc0ff0330000ffffffff0000f0e8d4aa552b170fc3c3a5a55a5a3c3c713369b24d96cc8e6666669966999999b296338e71cc694d5a5ac3c33c3ca5a555d417f00fe82baa000000ff00ffffffd4d4d4d4d4d4d4d45a5a5a5a5a5a5a5a6969696969696969000000000000000033333333333333333c3c3c3c3c3c3c3c17171717171717170000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[22]  = 4096'hb271d4e8172b8e4d3c3c3c3cc3c3c3c3f0556933cc96aa0f6666999966669999aaf0339669cc0f55a5a5a5a55a5a5a5a714d17d42be8b28e0000ffff0000ffff0f0ff0f0f0f00f0fc3a55a3c3c5aa5c38e8e717171718e8e99669966669966994d4db2b2b2b24d4d5ac33ca5a53cc35aaaaa55555555aaaaff00ff0000ff00ff8ee8b2d42b4d1771c33c3c3cc3c3c33caa33f069960fcc5599996699669966660f96aa33cc5569f05aa5a5a55a5a5aa5b2d47117e88e2b4dffff00ff00ff000033cccc3333cccc33a55a5aa5a55a5aa5d42b2bd4d42b2bd40066401910464411e81717e8e81717e8c33c3cc3c33c3cc369969669699696690000000000000000d44de8718e17b22b3cc33c3cc3c33cc3690f3355aaccf0969999996699666666335596f00f69aacca55aa5a55a5aa55a178ed44db22b71e8ffffff00ff000000f00ff00f0ff00ff05ac33ca5a53cc35a718e718e8e718e719999666666669999b24db24d4db24db23c5aa5c3c3a55a3c55aa55aaaa55aa55ffff00000000ffff17d48eb24d712be8c33cc33cc33cc33ccc69aaf00f559633669999669966669969330faa55f0cc965aa55aa55aa55aa52b17b2718e4de8d400ffff00ff0000ff3333333333333333a5a5a5a5a5a5a5a5d4d4d4d4d4d4d4d40000000000000000e8e8e8e8e8e8e8e8c3c3c3c3c3c3c3c369696969696969690000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[23]  = 4096'h175517f00fe8aae85aa5a53cc35a5aa5cc4dcc718e33b2336666666699999999697169b24d968e963cc3c3a55a3c3cc32b0f2b55aad4f0d400000000ffffffff698e719696718e693c5aa5c3c3a55a3c17aa55e8e855aa176666999999996666d40ff02b2bf00fd4a53cc35a5ac33ca5ccb24d33334db2cc0000ffffffff0000aaf01717e8e80f555a3c5aa55aa5c3a5b271cccc33338e4d99666666999999668eb2696996964d713ca53cc33cc35ac3f0552b2bd4d4aa0fff000000ffffff000ff0f00f0ff0f00f3cc3c33c3cc3c33c8e71718e8e71718e00664019104644114db2b24d4db2b24da55a5aa5a55a5aa5aa5555aaaa5555aa000000000000000017e8f055aa0f17e8a5a53ca55ac35a5acc33714db28ecc3366996666999966996996b2718e4d6996c3c3a5c33c5a3c3c2bd4550ff0aa2bd400ff0000ffff00ff7169968e8e966971a53cc35a5ac33ca55517e8aaaae817559966996666996699f0d42b0f0f2bd4f0c3a55a3c3c5aa5c34dcc33b2b233cc4dff00ff0000ff00ff0f17aa17e855e8f0c3a55a5aa5a55a3c8eccb2cc334d337199669966996699664d698e69967196b25ac33c3cc3c33ca5aa2bf02bd40fd455ff00ff00ff00ff000f0f0f0f0f0f0f0f3c3c3c3c3c3c3c3c8e8e8e8e8e8e8e8e00000000000000004d4d4d4d4d4d4d4da5a5a5a5a5a5a5a5aaaaaaaaaaaaaaaa0000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[24]  = 4096'h55cc17b24de833aac3c3c33cc33c3c3c4d2bccf00f33d4b29999666699996666711769aa5596e88e5a5a5aa55aa5a5a50f962b718ed469f0ffff0000ffff0000e8e817171717e8e8c3a55a3c3c5aa5c33333cccccccc3333996699666699669996966969696996965ac33ca5a53cc35ad4d42b2b2b2bd4d4ff00ff0000ff00ff33b25517e8aa4dcc3c3cc3c33c3cc3c3d4f04dcc33b20f2b6666996699669999e8aa7169968e5517a5a55a5aa5a55a5a69710f2bd4f08e960000ff00ff00ffffcc3333cccc3333cca55a5aa5a55a5aa52bd4d42b2bd4d42b006640191046441117e8e81717e8e817c33c3cc3c33c3cc39669699696696996000000000000000017aab2cc334d55e8c33c3cc33cc3c33cccb2f02bd40f4d336666669966999999698eaa17e85571965aa5a55aa55a5aa52bf07196698e0fd4000000ff00ffffff17e817e8e817e8175ac33ca5a53cc35acc33cc3333cc33cc999966666666999969966996966996693c5aa5c3c3a55a3c2bd42bd4d42bd42bffff00000000ffff4d173355aacce8b2c3c33cc33cc33c3c0fccd44db22b33f099666699669999665569e8718e1796aa5a5aa55aa55aa5a58e2b690ff096d471ff0000ff00ffff00cccccccccccccccca5a5a5a5a5a5a5a52b2b2b2b2b2b2b2b00000000000000001717171717171717c3c3c3c3c3c3c3c396969696969696960000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[25]  = 4096'hb271d4e8172b8e4dc35aa55aa55aa53cf0556933cc96aa0f9966996699669966aaf0339669cc0f555a3cc33cc33cc3a5714d17d42be8b28eff00ff00ff00ff00aad42b55552bd4aa3c5aa5c3c3a55a3cb269964d4d9669b266996699996699668e33cc7171cc338ea53cc35a5ac33ca5f017e80f0fe817f000ff00ffff00ff008ee8b2d42b4d1771a55ac3a55a3ca55aaa33f069960fcc5599669999666699660f96aa33cc5569f0c33c5ac33ca5c33cb2d47117e88e2b4dff00ffff0000ff00aa5555aaaa5555aa5aa5a55a5aa5a55ab24d4db2b24d4db200664019104644118e71718e8e71718e3cc3c33c3cc3c33cf00f0ff0f00f0ff00000000000000000d44de8718e17b22ba53c5a5aa5a5c35a690f3355aaccf0969966666699999966335596f00f69aaccc3a53c3cc3c35a3c178ed44db22b71e8ff000000ffffff002baa55d4d455aa2ba53cc35a5ac33ca596b24d69694db2966666999999996666cc8e713333718eccc3a55a3c3c5aa5c3e8f00f17170ff0e80000ffffffff000017d48eb24d712be8a5a5a5c33c5a5a5acc69aaf00f559633999999996666666669330faa55f0cc96c3c3c35aa53c3c3c2b17b2718e4de8d4ffffffff00000000aaaaaaaaaaaaaaaa5a5a5a5a5a5a5a5ab2b2b2b2b2b2b2b200000000000000008e8e8e8e8e8e8e8e3c3c3c3c3c3c3c3cf0f0f0f0f0f0f0f00000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[26]  = 4096'h0f0f2b17e8d4f0f03cc33c5aa5c33cc38e8e96cc3369717199666699669999664d4dcc699633b2b2a55aa53cc35aa55aaaaae82bd4175555ff0000ff00ffff004d55aab2b2aa554d3c5aa5c3c3a55a3c0f4db2f0f0b24d0f669966999966996655718eaaaa8e7155a53cc35a5ac33ca58e0ff07171f00f8e00ff00ffff00ff00f0170f2bd4f0e80f3c5a3c3cc3c3a5c371cc8e966971338e9999996699666666b2694dcc33b2964da53ca5a55a5ac35a552baae81755d4aaffffff00ff00000096696996966969965aa5a55a5aa5a55ae81717e8e81717e800664019104644112bd4d42b2bd4d42b3cc3c33c3cc3c33c33cccc3333cccc3300000000000000002bf0170ff0e80fd43cc35ac33ca53cc39671cc8e71338e696666996699669999ccb2694db2964d33a55a3c5aa5c3a55ae8552baa55d4aa170000ff00ff00ffffaa4db25555b24daaa53cc35a5ac33ca5b20ff04d4df00fb266669999999966668e55aa7171aa558ec3a55a3c3c5aa5c3f08e710f0f718ef00000ffffffff0000e82bf00ff00fd417a53c3c3cc3c3c35a3396718e718e69cc666699996666999996ccb24db24d3369c3a5a5a55a5a5a3cd4e855aa55aa172b0000ffff0000ffff96969696969696965a5a5a5a5a5a5a5ae8e8e8e8e8e8e8e800000000000000002b2b2b2b2b2b2b2b3c3c3c3c3c3c3c3c33333333333333330000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[27]  = 4096'h698eb233cc4d7196c3c3a5a55a5a3c3c17aaf0d42b0f55e86699999966666699d40faae81755f02b5a5ac3c33c3ca5a5ccb27169968e4d3300ffffff000000fff0cc330f0f33ccf0c3a55a3c3c5aa5c3712bd48e8ed42b716666999999996666b217e84d4de817b25ac33ca5a53cc35a559669aaaa6996550000ffffffff0000713369b24d96cc8e3ca5c3a55a3c5ac355d417f00fe82baa6699669966996699f0e8d4aa552b170fa5c35ac33ca53c5a4d69cc718e3396b200ff00ff00ff00ff718e8e71718e8e713cc3c33c3cc3c33c55aaaa5555aaaa550066401910464411f00f0ff0f00f0ff0a55a5aa5a55a5aa54db2b24d4db2b24d0000000000000000b296338e71cc694da53ca5c33c5ac35af0e8d4aa552b170f9999999966666666aa2be80ff017d455c3a5c35aa53c5a3c713369b24d96cc8effffffff0000000033f00fcccc0ff0335ac33ca5a53cc35ad4718e2b2b8e71d49966996666996699e8b24d17174db2e83c5aa5c3c3a55a3c6955aa9696aa5569ff00ff0000ff00ffccb27169968e4d335aa53cc33cc35aa52bf05517e8aa0fd4669966669999669917aaf0d42b0f55e83cc3a55aa55a3cc396714dcc33b28e6900ff0000ffff00ff71717171717171713c3c3c3c3c3c3c3c55555555555555550000000000000000f0f0f0f0f0f0f0f0a5a5a5a5a5a5a5a54d4d4d4d4d4d4d4d0000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[28]  = 4096'h55cc0f9669f033aaa55a5aa55aa5a55a4d2b8ee81771d4b2999999669966666671174d2bd4b2e88ec33c3cc33cc3c33c0f96aa33cc5569f0ffffff00ff00000096699669699669963c5aa5c3c3a55a3ce817e81717e817e899996666666699992bd42bd4d42bd42ba53cc35a5ac33ca533cc33cccc33cc33ffff00000000ffff3396550ff0aa69cca5a5a55aa55a5a5ad4e84d8e71b2172b6666999966669999e82b714db28ed417c3c3c33cc33c3c3c69330faa55f0cc960000ffff0000ffffe81717e8e81717e8c33c3cc3c33c3cc333cccc3333cccc33006640191046441196696996966969965aa5a55a5aa5a55ad42b2bd4d42b2bd400000000000000000faa96cc336955f05a5aa55aa55aa5a58eb2e82bd4174d7199666699669999664d8e2b17e8d471b23c3cc33cc33cc3c3aaf0339669cc0f55ff0000ff00ffff009696696969699696a53cc35a5ac33ca5e8e817171717e8e866996699996699662b2bd4d4d4d42b2bc3a55a3c3c5aa5c33333cccccccc333300ff00ffff00ff00690f3355aaccf0965a5aa5a55a5aa5a5178ed44db22b71e89999669966996666d44de8718e17b22b3c3cc3c33c3cc3c3ccaa690ff0965533ffff00ff00ff0000e8e8e8e8e8e8e8e8c3c3c3c3c3c3c3c33333333333333333000000000000000096969696969696965a5a5a5a5a5a5a5ad4d4d4d4d4d4d4d40000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[2]  = 4096'h144e0f9669f072283c81a542245a18c3488b8ee81771d1126600666699990099211d4d2bd4b2b884a542c318813c245a03a6aa33cc5565c000000000ffff00ff6115a88686a815613c5aa5c3c3a55a3c154c32a8a8324c156699669999669966d061860b0b8661d0a53cc35a5ac33ca54c0bd03232d00b4c00ff00ffff00ff007296140ff028694e18423ca55ac32481d1e8488e7112178b0066666699999900b82b214db284d41d2418a5c33c5a8142653303aa55c0cca600000000ffffff000db0b00d0db0b00d5aa5a55a5aa5a55a0e70700e0e70700e006640191046441145a2a24545a2a2453cc3c33c3cc3c33c8a51518a8a51518a00000000000000000f28964e726914f0a5c3428118243c5a8e12e88bd117487166996600009966994d842b1db8d421b2c35a18422481a53caac033a665cc035500ff000000ff00ffa8618615158661a8a53cc35a5ac33ca53215a84c4ca81532666699999999666686d00b61610bd086c3a55a3c3c5aa5c3d04c320b0b324cd00000ffffffff0000690f7214284ef09624a5183cc3815a42178ed148128b71e89966006699009966d44db821841db22b81c324a55a423c18ccaa6503c0a65533ff000000ff00ff000d0d0d0d0d0d0d0d5a5a5a5a5a5a5a5a0e0e0e0e0e0e0e0e000000000000000045454545454545453c3c3c3c3c3c3c3c8a8a8a8a8a8a8a8a0000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[29]  = 4096'h55b269d42b964daaa53cc33cc33cc35a4df0176996e80fb2996666996699996671aad433cc2b558ec3a55aa55aa55a3c0f71cc17e8338ef0ff0000ff00ffff00aa2bd45555d42baa5ac33ca5a53cc35ab296694d4d6996b266669999999966668ecc33717133cc8e3c5aa5c3c3a55a3cf0e8170f0f17e8f00000ffffffff00004dd4556996aa2bb2c33ca5c33c5ac33c0f694d17e8b296f09999996699666666553371d42b8eccaa5aa5c35aa53c5aa58e170fcc33f0e871ffffff00ff00000096696996966969963cc3c33c3cc3c33ce81717e8e81717e8ff99bfe6efb9bbee2bd4d42b2bd4d42ba55a5aa5a55a5aa533cccc3333cccc33ffffffffffffffff69aad4b24d2b5596c35a3c3cc3c3a53c17b269f00f964de86666996699669999d48e33aa55cc712b5a3ca5a55a5ac3a5ccf017718ee80f330000ff00ff00ffffd4aa552b2b55aad43c5aa5c3c3a55a3c69b24d96964db2699966996666996699338e71cccc718e33a53cc35a5ac33ca517f00fe8e80ff017ff00ff0000ff00ff2b694d55aab296d4c3c3c3a55a3c3c3c96170f4db2f0e8696666999966669999ccd455718eaa2b335a5a5ac33ca5a5a5e8cc8e0ff07133170000ffff0000ffff96969696969696963c3c3c3c3c3c3c3ce8e8e8e8e8e8e8e8ffffffffffffffff2b2b2b2b2b2b2b2ba5a5a5a5a5a5a5a53333333333333333ffffffffffffffff >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[30]  = 4096'h55aa55aa55aa55aaa55aa55aa55aa55a4db24db24db24db29966996699669966718e718e718e718ec33cc33cc33cc33c0ff00ff00ff00ff0ff00ff00ff00ff00aa55aa5555aa55aa5aa55aa5a55aa55ab24db24d4db24db266996699996699668e718e71718e718e3cc33cc3c33cc33cf00ff00f0ff00ff000ff00ffff00ff0055aa5555aaaa55aaa55aa5a55a5aa55a4db24d4db2b24db29966999966669966718e71718e8e718ec33cc3c33c3cc33c0ff00f0ff0f00ff0ff00ffff0000ff00aa5555aaaa5555aa5aa5a55a5aa5a55ab24d4db2b24d4db266999966669999668e71718e8e71718e3cc3c33c3cc3c33cf00f0ff0f00f0ff000ffff0000ffff0055aaaaaa555555aaa55a5a5aa5a5a55a4db2b2b24d4d4db29966666699999966718e8e8e7171718ec33c3c3cc3c3c33c0ff0f0f00f0f0ff0ff000000ffffff00aaaa55555555aaaa5a5aa5a5a5a55a5ab2b24d4d4d4db2b266669999999966668e8e717171718e8e3c3cc3c3c3c33c3cf0f00f0f0f0ff0f00000ffffffff000055555555aaaaaaaaa5a5a5a55a5a5a5a4d4d4d4db2b2b2b29999999966666666717171718e8e8e8ec3c3c3c33c3c3c3c0f0f0f0ff0f0f0f0ffffffff00000000aaaaaaaaaaaaaaaa5a5a5a5a5a5a5a5ab2b2b2b2b2b2b2b266666666666666668e8e8e8e8e8e8e8e3c3c3c3c3c3c3c3cf0f0f0f0f0f0f0f00000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[31]  = 4096'h55aa55aa55aa55aaa55aa55aa55aa55a4db24db24db24db29966996699669966718e718e718e718ec33cc33cc33cc33c0ff00ff00ff00ff0ff00ff00ff00ff00aa55aa5555aa55aa5aa55aa5a55aa55ab24db24d4db24db266996699996699668e718e71718e718e3cc33cc3c33cc33cf00ff00f0ff00ff000ff00ffff00ff0055aa5555aaaa55aaa55aa5a55a5aa55a4db24d4db2b24db29966999966669966718e71718e8e718ec33cc3c33c3cc33c0ff00f0ff0f00ff0ff00ffff0000ff00aa5555aaaa5555aa5aa5a55a5aa5a55ab24d4db2b24d4db266999966669999668e71718e8e71718e3cc3c33c3cc3c33cf00f0ff0f00f0ff000ffff0000ffff0055aaaaaa555555aaa55a5a5aa5a5a55a4db2b2b24d4d4db29966666699999966718e8e8e7171718ec33c3c3cc3c3c33c0ff0f0f00f0f0ff0ff000000ffffff00aaaa55555555aaaa5a5aa5a5a5a55a5ab2b24d4d4d4db2b266669999999966668e8e717171718e8e3c3cc3c3c3c33c3cf0f00f0f0f0ff0f00000ffffffff000055555555aaaaaaaaa5a5a5a55a5a5a5a4d4d4d4db2b2b2b29999999966666666717171718e8e8e8ec3c3c3c33c3c3c3c0f0f0f0ff0f0f0f0ffffffff00000000aaaaaaaaaaaaaaaa5a5a5a5a5a5a5a5ab2b2b2b2b2b2b2b266666666666666668e8e8e8e8e8e8e8e3c3c3c3c3c3c3c3cf0f0f0f0f0f0f0f00000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[3]  = 4096'h56d517f00fe8ab6ac37e5abddba5e73cc96dcc718e33b6936600666699990099397369b24d96ce9c5abd3ce77ec3dba5271f2b55aad4f8e400000000ffff00ff9eea57797957ea9ec3a55a3c3c5aa5c3eab3cd5757cdb3ea99996666666699992f9e79f4f4799e2f5ac33ca5a53cc35ab3f42fcdcd2ff4b3ffff00000000ffffabf05617e86a0fd5e7bdc35aa53cdb7eb671c9cc33938e6d0066666699999900ceb23969969c4d73dbe75a3cc3a57ebdf855272bd4e4aa1f00000000ffffff000db0b00d0db0b00dc33c3cc3c33c3cc30e70700e0e70700e006640191046441145a2a24545a2a2455aa5a55a5aa5a55a8a51518a8a51518a0000000000000000176af0d5ab0f56e85a3cbd7ee7dbc3a5cc93716db68ec9336699660000996699699cb273ce4d39963ca5e7bddb7e5ac32be4551ff8aa27d400ff000000ff00ff579e79eaea799e575ac33ca5a53cc35acdea57b3b357eacd6699669999669966792ff49e9ef42f793c5aa5c3c3a55a3c2fb3cdf4f4cdb32f00ff00ffff00ff000f17ab566ad5e8f0db5ae7c33c7ea5bd8eccb6c9936d337199660066990099664d69ce399c7396b27e3cdb5aa5bdc3e7aa2bf827e41fd455ff000000ff00ff000d0d0d0d0d0d0d0dc3c3c3c3c3c3c3c30e0e0e0e0e0e0e0e000000000000000045454545454545455a5a5a5a5a5a5a5a8a8a8a8a8a8a8a8a0000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[4]  = 4096'he9b2e84db2174d97a5a53c3cc3c35a5a37f0330ff0cc0fec9900666699990066d6aa9655aa69556bc3c3a5a55a5a3c3cdc71d48e712b8e3bff000000ffff0000cc0ff03333f00fccc3a55a3c3c5aa5c32b8e71d4d4718e2b9999666666669999174db2e8e8b24d175ac33ca5a53cc35a96aa55696955aa96ffff00000000ffff4d4de9e81797b2b25a3ca53cc35ac3a50f0f3733ccecf0f000669966996699005555d696696baaaa3ca5c3a55a3c5ac38e8edcd42b3b71710000ff00ff00ff008c31318c8c31318cc33c3cc3c33c3cc32a54542a2a54542a006640191046441107e0e00707e0e0075aa5a55a5aa5a55a92494992924949920000000000000000e8974db24db2e9173c5a3ca55ac3a5c333ec0ff00ff037cc6666660000999999966b55aa55aad669a53ca5c33c5ac35ad43b8e718e71dc2b0000000000fffffff0cc330f0f33ccf05ac33ca5a53cc35a712bd48e8ed42b716699669999669966b217e84d4de817b23c5aa5c3c3a55a3c559669aaaa69965500ff00ffff00ff00b2e84de997b2174dc33c5aa55aa5c33cf0330f37ecf0cc0f9966009966009966aa9655d66baa69555aa53cc33cc35aa571d48edc3b712b8eff0000ff0000ff008c8c8c8c8c8c8c8cc3c3c3c3c3c3c3c32a2a2a2a2a2a2a2a000000000000000007070707070707075a5a5a5a5a5a5a5a92929292929292920000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[5]  = 4096'h96714dcc33b28e695a3ca53cc35ac3a5e8550f2bd4f0aa1766009999666600992bf05517e8aa0fd43ca5c3a55a3c5ac3334d8e966971b2cc0000ffff000000ff2b8e71d4d4718e2b3c5aa5c3c3a55a3c96aa55696955aa969966996666996699cc0ff03333f00fcca53cc35a5ac33ca5e8b24d17174db2e8ff00ff0000ff00ff8ecc964db2693371c33c5aa55aa5c33caa2be80ff017d45500996699669966000f172b55aad4e8f05aa53cc33cc35aa5b296338e71cc694d00ff00ff00ff0000318c8c31318c8c31a55a5aa5a55a5aa5542a2a54542a2a540066401910464411e00707e0e00707e0c33c3cc3c33c3cc3499292494992924900000000000000004d69cc718e3396b2a5a53c3cc3c35a5a0f172b55aad4e8f0999999000066666655d417f00fe82baac3c3a5a55a5a3c3c8ecc964db2693371ffffff0000000000712bd48e8ed42b71a53cc35a5ac33ca5559669aaaa6996559999666666669999f0cc330f0f33ccf0c3a55a3c3c5aa5c34de817b2b217e84dffff00000000ffff334d8e966971b2ccc3a5c35aa53c5a3cd40faae81755f02b6699006699006699e8550f2bd4f0aa175ac35a3cc3a53ca5698eb233cc4d719600ff0000ff0000ff3131313131313131a5a5a5a5a5a5a5a554545454545454540000000000000000e0e0e0e0e0e0e0e0c3c3c3c3c3c3c3c349494949494949490000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[6]  = 4096'h0f0f3333ccccf0f0a55aa5a55a5aa55a8e8ed4d42b2b717166ff66996699ff994d4de8e81717b2b2c33cc3c33c3cc33caaaa69699696555500ff00ff00ffffffb24db24d4db24db23c5aa5c3c3a55a3cf00ff00f0ff00ff09999666666669999aa55aa5555aa55aaa53cc35a5ac33ca5718e718e8e718e71ffff00000000fffff0330f33ccf0cc0fa5a5a5a55a5a5a5a71d48ed42b712b8eff996666999966ffb2e84de817b2174dc3c3c3c33c3c3c3c5569aa69965596aaffff0000ffff00ff57eaea5757eaea57c33c3cc3c33c3cc3cdb3b3cdcdb3b3cd0066401910464411799e9e79799e9e795aa5a55a5aa5a55a2ff4f42f2ff4f42f000000000000000033f0330ff0cc0fcca55aa55aa55aa55ad471d48e712b8e2b669999ffff666699e8b2e84db2174d17c33cc33cc33cc33c695569aa5596aa9600ffffffff0000ffb2b24d4d4d4db2b2a53cc35a5ac33ca5f0f00f0f0f0ff0f06699669999669966aaaa55555555aaaac3a55a3c3c5aa5c371718e8e8e8e717100ff00ffff00ff00cc33f00ff00fcc335aa5a5a55a5a5aa52bd4718e718e2bd46666ff6699ff999917e8b24db24d17e83cc3c3c33c3c3cc3966955aa55aa96690000ff00ffffffff5757575757575757c3c3c3c3c3c3c3c3cdcdcdcdcdcdcdcd000000000000000079797979797979795a5a5a5a5a5a5a5a2f2f2f2f2f2f2f2f0000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[7]  = 4096'h55d417f00fe82baa3cc3a53cc35a3cc34d69cc718e3396b26699666699996699713369b24d96cc8ea55ac3a55a3ca55a0f172b55aad4e8f000ff0000ffff00ff694db29696b24d69c3a55a3c3c5aa5c3170ff0e8e8f00f176699669999669966d455aa2b2baa55d45ac33ca5a53cc35acc8e713333718ecc00ff00ffff00ff002bf05517e8aa0fd43c3c3ca55ac3c3c396714dcc33b28e696666666699999999ccb27169968e4d33a5a5a5c33c5a5a5ae8550f2bd4f0aa1700000000ffffffff4db2b24d4db2b24d5aa5a55a5aa5a55a0ff0f00f0ff0f00f006640191046441155aaaa5555aaaa553cc3c33c3cc3c33c8e71718e8e71718e000000000000000017aaf0d42b0f55e8a5c33cc33cc33c5accb27169968e4d336699669966996699698eb233cc4d7196c35aa55aa55aa53c2bf05517e8aa0fd400ff00ff00ff00ffb269964d4d9669b25ac33ca5a53cc35af017e80f0fe817f06666999999996666aad42b55552bd4aa3c5aa5c3c3a55a3c71cc338e8e33cc710000ffffffff00000f172b55aad4e8f0c3a53c3cc3c35a3c8ecc964db269337199666666999999664d69cc718e3396b25ac3a5a55a5a3ca5aa2be80ff017d455ff000000ffffff004d4d4d4d4d4d4d4d5a5a5a5a5a5a5a5a0f0f0f0f0f0f0f0f000000000000000055555555555555553c3c3c3c3c3c3c3c8e8e8e8e8e8e8e8e0000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
  assign \coef[8]  = 4096'h2b17b2718e4de8d45a5ac3a55a3ca5a596ccf055aa0f33699999669966996666cc69aaf00f5596333c3c5ac33ca5c3c3e82b714db28ed417ffff00ff00ff000033e817cccc17e833c3a55a3c3c5aa5c3d433cc2b2bcc33d46666999999996666e8966917176996e85ac33ca5a53cc35a69d42b96962bd4690000ffffffff0000e8712bb24dd48e17a5a55ac33ca55a5a335596f00f69aacc669999669966669996f0ccaa55330f69c3c33c5aa5c33c3cd44de8718e17b22b00ffff00ff0000ffd42b2bd4d42b2bd43cc3c33c3cc3c33c6996966969969669006640191046441133cccc3333cccc33a55a5aa5a55a5aa517e8e81717e8e8170000000000000000b2d47117e88e2b4dc3a5a55aa55a5a3cf06955cc33aa960f6666999966669999aa33f069960fcc555ac3c33cc33c3ca571174d2bd4b2e88e0000ffff0000ffff1733cce8e8cc33175ac33ca5a53cc35accd42b33332bd4cc996699666699669969e817969617e8693c5aa5c3c3a55a3c2b6996d4d496692bff00ff0000ff00ff8eb2e82bd4174d715ac3a55aa55a3ca5aaf0339669cc0f5566666699669999990faa96cc336955f03c5ac33cc33ca5c3b271d4e8172b8e4d000000ff00ffffffd4d4d4d4d4d4d4d43c3c3c3c3c3c3c3c696969696969696900000000000000003333333333333333a5a5a5a5a5a5a5a517171717171717170000000000000000 >> { \v[2] , \v[1] , \v[0] , \u[2] , \u[1] , \u[0] , \y[2] , \y[1] , \y[0] , \x[2] , \x[1] , \x[0]  };
endmodule
