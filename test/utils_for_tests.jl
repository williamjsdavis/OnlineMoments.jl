# Test utility functions

"Load some small test data from a local file"
function load_small_data()
    filename = "./X_data_small.jld2"
    raw_data = load(filename)
    return raw_data["X"]
end

"Transforms a static array to a streaming-like object"
function stream_data(data)
    data_copy = deepcopy(data)
    data_stream = () -> popfirst!(data_copy)
    return data_stream
end