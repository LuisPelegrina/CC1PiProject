//
// Created by Luis Pelegrina Gutiérrez on 30/6/24.
//
#include "../Includes.h"

map<int, int> resonance_index_map = {
        {0 ,0},
        {202228 ,50},
        {202218 ,49},
        {202118 ,48},
        {201118 ,47},

        {222224 ,46},
        {222214 ,45},
        {222114 ,44},
        {221114 ,43},

        {212226 ,42},
        {212216 ,41},
        {212116 ,40},
        {211116 ,39},

        {212212 ,38},
        {212112 ,37},

        {102216 ,36},
        {102116 ,35},

        {222222 ,34},
        {222212 ,33},
        {222112 ,32},
        {221112 ,31},

        {112214 ,30},
        {112114 ,29},

        {132212 ,28},
        {132112 ,27},

        {122224 ,26},
        {122214 ,25},
        {122114 ,24},
        {121114 ,23},

        {102212 ,22},
        {102112 ,21},

        {112222 ,20},
        {112212 ,19},
        {112112 ,18},
        {111112 ,17},

        {202114 ,16},
        {202214 ,15},

        {202212 ,14},
        {202112 ,13},

        {202216 ,12},
        {202116 ,11},

        {102214 ,10},
        {102114 ,9},

        {212224 ,8},
        {212214 ,7},
        {212114 ,6},
        {211114 ,5},

        {2224 ,4},
        {2214 ,3},
        {2114 ,2},
        {1114 ,1}
};

map<int, string> resonance_type_delta_or_N_map = {
        {0 ,"unclasified"},
        {50, "delta"},
        {49, "delta"},
        {48, "delta"},
        {47, "delta"},

        {46, "delta"},
        {45, "delta"},
        {44, "delta"},
        {43, "delta"},

        {42, "delta"},
        {41, "delta"},
        {40, "delta"},
        {39, "delta"},

        {38, "np"},
        {37, "np"},

        {36, "np"},
        {35, "np"},

        {34, "delta"},
        {33, "delta"},
        {32, "delta"},
        {31, "delta"},

        {30, "np"},
        {29, "np"},

        {28, "np"},
        {27, "np"},

        {26, "delta"},
        {25, "delta"},
        {24, "delta"},
        {23, "delta"},

        {22, "np"},
        {21, "np"},

        {20, "delta"},
        {19, "delta"},
        {18, "delta"},
        {17, "delta"},

        {16, "np"},
        {15, "np"},

        {14, "np"},
        {13, "np"},

        {12, "np"},
        {11, "np"},

        {10, "np"},
        {9, "np"},

        {8, "delta"},
        {7, "delta"},
        {6, "delta"},
        {5, "delta"},

        {4, "delta"},
        {3, "delta"},
        {2, "delta"},
        {1, "delta"}
};


map<int, string> resonance_type_name_map = {
        {0 ,"unclasified"},
        {50, "Delta (1950)"},
        {49, "Delta (1950)"},
        {48, "Delta (1950)"},
        {47, "Delta (1950)"},

        {46, "Delta (1920)"},
        {45, "Delta (1920)"},
        {44, "Delta (1920)"},
        {43, "Delta (1920)"},

        {42, "Delta (1905)"},
        {41, "Delta (1905)"},
        {40, "Delta (1905)"},
        {39, "Delta (1905)"},

        {38, "N (1710)"},
        {37, "N (1710)"},

        {36, "N (1675)"},
        {35, "N (1675)"},

        {34, "Delta (1910)"},
        {33, "Delta (1910)"},
        {32, "Delta (1910)"},
        {31, "Delta (1910)"},

        {30, "N (1700)"},
        {29, "N (1700)"},

        {28, "N (1650)"},
        {27, "N (1650)"},

        {26, "Delta (1700)"},
        {25, "Delta (1700)"},
        {24, "Delta (1700)"},
        {23, "Delta (1700)"},

        {22, "N (1535)"},
        {21, "N (1535)"},

        {20, "Delta (1620)"},
        {19, "Delta (1620)"},
        {18, "Delta (1620)"},
        {17, "Delta (1620)"},

        {16, "N (1720)"},
        {15, "N (1720)"},

        {14, "N (1440)"},
        {13, "N (1440)"},

        {12, "N (1680)"},
        {11, "N (1680)"},

        {10, "N (1520)"},
        {9,  "N (1520)"},

        {8, "Delta (1600)"},
        {7, "Delta (1600)"},
        {6, "Delta (1600)"},
        {5, "Delta (1600)"},

        {4, "Delta (1232)"},
        {3, "Delta (1232)"},
        {2, "Delta (1232)"},
        {1, "Delta (1232)"}
};

string GetTitleRes(int hIndex) {
    string sTitle;
    switch (hIndex) {
        case 50:
            sTitle = "#Delta^{++}(1950) resonance; E_{Nu}; Events";
            break;
        case 49:
            sTitle = "#Delta^{+}(1950) resonance; E_{Nu}; Events";
            break;
        case 48:
            sTitle = "#Delta^{0}(1950) resonance; E_{Nu}; Events";
            break;
        case 47:
            sTitle = "#Delta^{-}(1950) resonance; E_{Nu}; Events";
            break;

        case 46:
            sTitle = "#Delta^{++}(1920) resonance; E_{Nu}; Events";
            break;
        case 45:
            sTitle = "#Delta^{+}(1920) resonance; E_{Nu}; Events";
            break;
        case 44:
            sTitle = "#Delta^{0}(1920) resonance; E_{Nu}; Events";
            break;
        case 43:
            sTitle = "#Delta^{-}(1920) resonance; E_{Nu}; Events";
            break;

        case 42:
            sTitle = "#Delta^{++}(1905) resonance; E_{Nu}; Events";
            break;
        case 41:
            sTitle = "#Delta^{+}(1905) resonance; E_{Nu}; Events";
            break;
        case 40:
            sTitle = "#Delta^{0}(1905) resonance; E_{Nu}; Events";
            break;
        case 39:
            sTitle = "#Delta^{-}(1905) resonance; E_{Nu}; Events";
            break;

        case 38:
            sTitle = "p(1710) resonance; E_{Nu}; Events";
            break;
        case 37:
            sTitle = "n(1710) resonance; E_{Nu}; Events";
            break;

        case 36:
            sTitle = "p(1675) resonance; E_{Nu}; Events";
            break;
        case 35:
            sTitle = "n(1675) resonance; E_{Nu}; Events";
            break;

        case 34:
            sTitle = "#Delta^{++}(1910) resonance; E_{Nu}; Events";
            break;
        case 33:
            sTitle = "#Delta^{+}(1910) resonance; E_{Nu}; Events";
            break;
        case 32:
            sTitle = "#Delta^{0}(1910) resonance; E_{Nu}; Events";
            break;
        case 31:
            sTitle = "#Delta^{-}(1910) resonance; E_{Nu}; Events";
            break;

        case 30:
            sTitle = "p(1700) resonance; E_{Nu}; Events";
            break;
        case 29:
            sTitle = "n(1700) resonance; E_{Nu}; Events";
            break;

        case 28:
            sTitle = "p(1650) resonance; E_{Nu}; Events";
            break;
        case 27:
            sTitle = "n(1650) resonance; E_{Nu}; Events";
            break;

        case 26:
            sTitle = "#Delta^{++}(1700) resonance; E_{Nu}; Events";
            break;
        case 25:
            sTitle = "#Delta^{+}(1700) resonance; E_{Nu}; Events";
            break;
        case 24:
            sTitle = "#Delta^{0}(1700) resonance; E_{Nu}; Events";
            break;
        case 23:
            sTitle = "#Delta^{-}(1700) resonance; E_{Nu}; Events";
            break;

        case 22:
            sTitle = "p(1535) resonance; E_{Nu}; Events";
            break;
        case 21:
            sTitle = "n(1535) resonance; E_{Nu}; Events";
            break;

        case 20:
            sTitle = "#Delta^{++}(1620) resonance; E_{Nu}; Events";
            break;
        case 19:
            sTitle = "#Delta^{+}(1620) resonance; E_{Nu}; Events";
            break;
        case 18:
            sTitle = "#Delta^{0}(1620) resonance; E_{Nu}; Events";
            break;
        case 17:
            sTitle = "#Delta^{-}(1620) resonance; E_{Nu}; Events";
            break;

        case 16:
            sTitle = "p(1720) resonance; E_{Nu}; Events";
            break;
        case 15:
            sTitle = "n(1720) resonance; E_{Nu}; Events";
            break;

        case 14:
            sTitle = "p(1440) resonance; E_{Nu}; Events";
            break;
        case 13:
            sTitle = "n(1440) resonance; E_{Nu}; Events";
            break;

        case 12:
            sTitle = "p(1680) resonance; E_{Nu}; Events";
            break;
        case 11:
            sTitle = "n(1680) resonance; E_{Nu}; Events";
            break;

        case 10:
            sTitle = "p(1520) resonance; E_{Nu}; Events";
            break;
        case 9:
            sTitle = "n(1520) resonance; E_{Nu}; Events";
            break;

        case 8:
            sTitle = "#Delta^{++}(1600) resonance; E_{Nu}; Events";
            break;
        case 7:
            sTitle = "#Delta^{+}(1600) resonance; E_{Nu}; Events";
            break;
        case 6:
            sTitle = "#Delta^{0}(1600) resonance; E_{Nu}; Events";
            break;
        case 5:
            sTitle = "#Delta^{-}(1600) resonance; E_{Nu}; Events";
            break;

        case 4:
            sTitle = "#Delta^{++}(1232) resonance; E_{Nu}; Events";
            break;
        case 3:
            sTitle = "#Delta^{+}(1232) resonance; E_{Nu}; Events";
            break;
        case 2:
            sTitle = "#Delta^{0}(1232) resonance; E_{Nu}; Events";
            break;
        case 1:
            sTitle = "#Delta^{-}(1232) resonance; E_{Nu}; Events";
            break;

        case 0:
            sTitle = "No Resonance; E_{Nu}; Events";
            break;

    }
    return sTitle;
}

string GetTitletog(int hIndex) {
    string sTitle;
    switch (hIndex) {
        case 0:
            sTitle = " Unknown resonance; E_{#nu} [GeV]; Events";
            break;
        case 1:
            sTitle = "#Delta resonance; E_{#nu} [GeV]; Events";
            break;
        case 2:
            sTitle = "N resonance; E_{#nu} [GeV]; Events";
            break;
    }
    return sTitle;
}
