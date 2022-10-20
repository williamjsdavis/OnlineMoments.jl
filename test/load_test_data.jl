# Loading test data from files

function load_small_data()
    filename = "./X_data_small.jld2"
    raw_data = load(filename)
    return raw_data["X"]
end
