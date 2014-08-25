/*
 * NPG-explorer, Nucleotide PanGenome explorer
 * Copyright (C) 2014 Boris Nagaev
 *
 * See the LICENSE file for terms of use.
 */

#define NPGE_SCRIPT(...) #__VA_ARGS__

const char* meta_lua = NPGE_SCRIPT(

function terminal()
    while true do
        io.write("> ")
        local line = io.read()
        if not line then
            break
        end
        local chunk, load_error = loadstring(line)
        if chunk then
            local success, pcall_result = pcall(chunk)
            print(pcall_result)
        else
            print(load_error)
        end;
    end
end

go = function(x) return meta:get_opt(x) end

function new_p(name)
    return meta:get_plain(name)
end

function run(name, opts)
    local p = new_p(name)
    opts = opts or ""
    p:set_options(opts, meta:placeholder_processor())
    p:run()
    Processor.delete(p)
end

function run_main(name)
    local args = {}
    for word in main_args:gmatch("%S+") do
        table.insert(args, word)
    end
    // remove first (program name)
    table.remove(args, 1)
    local args_s = ''
    for i, v in next, args do
        args_s = args_s..' '..v
    end
    run(name, args_s);
end

function block_set()
    return meta:placeholder_processor():block_set()
end

function other()
    return meta:placeholder_processor():other()
end

function register_p(name, returner)
    local returner1 = function()
        local p = returner()
        p:set_key(name)
        return p
    end
    meta:set_returner(returner1, name)
end

);

