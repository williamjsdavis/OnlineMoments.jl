# Test utility functions

"Transforms a static array to a streaming-like object"
function stream_data(data)
    data_copy = deepcopy(data)
    data_stream = () -> popfirst!(data_copy)
    return data_stream
end