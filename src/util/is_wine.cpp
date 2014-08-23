/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#if defined(_WIN32) || defined(__WIN32__)

#include <windows.h>

namespace npge {

bool is_wine() {
    // http://pgolub.wordpress.com/2009/05/18/on-wine-or-not-on-wine/
    HINSTANCE i = LoadLibrary("ntdll.dll");
    if (GetProcAddress(i, "wine_get_version")) {
        return true;
    } else {
        return false;
    }
}

}

#else

namespace npge {

bool is_wine() {
    return false;
}

}

#endif

