//
// Created by Luis Pelegrina Gutiérrez on 30/6/24.
//
#include "../../Includes.h"

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
