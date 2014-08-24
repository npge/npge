// http://stackoverflow.com/a/5841066

#define NPGE_SCRIPT(...) #__VA_ARGS__

const char* terminal_lua = NPGE_SCRIPT(

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

);

