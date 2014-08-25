/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#define NPGE_SCRIPT(...) #__VA_ARGS__

const char* meta_lua = NPGE_SCRIPT(

go = function(x) return meta:get_opt(x) end

function new_p(name)
    return meta:get_plain(name)
end

function register_p(name, returner)
    meta:set_returner(returner, name)
end

);

