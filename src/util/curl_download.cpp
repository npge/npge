/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <curl/curl.h>
#include <curl/easy.h>

#include "curl_download.hpp"
#include "name_to_stream.hpp"

namespace npge {

size_t write_callback(char* ptr, size_t size,
                      size_t nmemb, void* o0) {
    std::ostream* o = reinterpret_cast<std::ostream*>(o0);
    o->write(ptr, size * nmemb);
    return size * nmemb;
}

bool download_file(const std::string& url,
                   const std::string& out_fname) {
    typedef boost::shared_ptr<std::ostream> OPtr;
    OPtr out = name_to_ostream(out_fname);
    typedef boost::shared_ptr<CURL> CurlPtr;
    CurlPtr curl(curl_easy_init(), curl_easy_cleanup);
    curl_easy_setopt(curl.get(), CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl.get(), CURLOPT_WRITEFUNCTION,
                     write_callback);
    curl_easy_setopt(curl.get(), CURLOPT_WRITEDATA, out.get());
    curl_easy_setopt(curl.get(), CURLOPT_FOLLOWLOCATION, 1);
    return (curl_easy_perform(curl.get()) == CURLE_OK);
}

}

