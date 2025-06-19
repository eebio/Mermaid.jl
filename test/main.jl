@testitem "fullname" begin
    @test Mermaid.fullname(ConnectedVariable("comp", "var", 1:5, [1, 3, 5, 7])) == "comp[[1, 3, 5, 7]].var[1:5]"
    @test Mermaid.fullname(ConnectedVariable("cmp", "var", 1:5, nothing)) == "cmp.var[1:5]"
    @test Mermaid.fullname(ConnectedVariable("cp", "vr", nothing, [1, 3, 5, 7])) == "cp[[1, 3, 5, 7]].vr"
    @test Mermaid.fullname(ConnectedVariable("comp", "var", nothing, nothing)) == "comp.var"
end
