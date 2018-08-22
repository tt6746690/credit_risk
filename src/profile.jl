using Coverage

infos = analyze_malloc("src")
for info in infos
    println(info)
end
