/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#include <istream>
#include <ostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "download_file.hpp"
#include "name_to_stream.hpp"
#include "throw_assert.hpp"
#include "cast.hpp"
#include "Exception.hpp"

namespace npge {

bool download_file(const std::string& url,
                   const std::string& out_fname) {
    using namespace boost::algorithm;
    std::string http("http://");
    ASSERT_TRUE(starts_with(url, http));
    size_t slash_pos = url.find('/', http.size());
    ASSERT_TRUE(slash_pos != std::string::npos);
    int server_size = slash_pos - http.size();
    std::string server = url.substr(http.size(), server_size);
    std::string path = url.substr(slash_pos);
    // based on boost_asio/example/http/client/sync_client.cpp
    using boost::asio::ip::tcp;
    typedef boost::shared_ptr<std::ostream> OPtr;
    OPtr out = name_to_ostream(out_fname);
    std::ostream& o = *out;
    //
    boost::asio::io_service io_service;
    // Get a list of endpoints corresponding
    // to the server name.
    tcp::resolver resolver(io_service);
    tcp::resolver::query query(server, "http");
    tcp::resolver::iterator endpoint_iterator =
        resolver.resolve(query);
    // Try each endpoint until we
    // successfully establish a connection.
    tcp::socket socket(io_service);
    boost::asio::connect(socket, endpoint_iterator);
    // Form the request.
    // We specify the "Connection: close" header so that the
    // server will close the socket after transmitting
    // the response. This will
    // allow us to treat all data up until the EOF
    // as the content.
    boost::asio::streambuf request;
    std::ostream request_stream(&request);
    request_stream << "GET " << path << " HTTP/1.0\r\n";
    request_stream << "Host: " << server << "\r\n";
    request_stream << "Accept: */*\r\n";
    request_stream << "Connection: close\r\n\r\n";
    // Send the request.
    boost::asio::write(socket, request);
    // Read the response status line.
    // The response streambuf will automatically
    // grow to accommodate the entire line.
    // The growth may be limited by passing
    // a maximum size to the streambuf constructor.
    boost::asio::streambuf response;
    boost::asio::read_until(socket, response, "\r\n");
    // Check that response is OK.
    std::istream response_stream(&response);
    std::string http_version;
    response_stream >> http_version;
    unsigned int status_code;
    response_stream >> status_code;
    std::string status_message;
    std::getline(response_stream, status_message);
    if (!response_stream ||
            http_version.substr(0, 5) != "HTTP/") {
        return false;
    }
    if (status_code != 200) {
        return false;
    }
    // Read the response headers,
    // which are terminated by a blank line.
    boost::asio::read_until(socket, response, "\r\n\r\n");
    // Process the response headers.
    std::string header;
    while (std::getline(response_stream, header) &&
            header != "\r") {
    }
    // Write whatever content we already have to output.
    if (response.size() > 0) {
        o << &response;
    }
    // Read until EOF, writing data to output as we go.
    boost::system::error_code error;
    while (boost::asio::read(
                socket, response,
                boost::asio::transfer_at_least(1), error)) {
        o << &response;
    }
    if (error != boost::asio::error::eof) {
        return false;
    }
    return true;
}

}

