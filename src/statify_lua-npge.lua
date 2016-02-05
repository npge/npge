--
-- NPG-explorer, Nucleotide PanGenome explorer
-- Copyright (C) 2012-2016 Boris Nagaev
--
-- See the LICENSE file for terms of use.
--

local npge_dir = assert(arg[1]) .. '/'
local lua_npge_dir = npge_dir .. 'src/lua-npge/'

local function readFile(path)
    local f = io.open(path)
    local content = f:read('*all')
    f:close()
    return content
end

dofile(lua_npge_dir .. 'lua-npge-dev-1.rockspec')
for m, file in pairs(_G.build.modules) do
    if m ~= 'npge.cpp' then
        assert(type(file) == 'string')
        local content = readFile(lua_npge_dir .. file)
        print(([[package.preload[%q] = function(...)
            %s
        end]]):format(m, content))
    end
end

print("npge = require 'npge'")
