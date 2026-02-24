@testitem "state control" begin
    function f(x)
        return x^2
    end
    comp = TimeIndependentComponent("squared", f, 2.0)
    int = init(comp)

    @test getstate(int, ConnectedVariable("squared.#state")) == 4.0
    @test getstate(int) == 4.0
    setstate!(int, ConnectedVariable("squared.#state"), 3.0)
    @test getstate(int, ConnectedVariable("squared.#state")) == 9.0
    @test getstate(int) == 9.0
    setstate!(int, 4.0)
    setstate!(int, ConnectedVariable("squared.#state"), 5.0)
    @test getstate(int, ConnectedVariable("squared.#state")) == 25.0
    @test getstate(int) == 25.0
end

@testitem "efficiency" begin
    # Makes sure we don't call the function more than necessary.
    function f(x)
        sleep(0.1)
        return x^2
    end
    comp = TimeIndependentComponent("squared", f, 2.0)
    int = init(comp)

    function timer(int)
        setstate!(int, ConnectedVariable("squared.#state"), 3.0)
        setstate!(int, ConnectedVariable("squared.#state"), 4.0)
        setstate!(int, ConnectedVariable("squared.#state"), 5.0)
        getstate(int, ConnectedVariable("squared.#state")) # required computation
        getstate(int) # should not require computation
        getstate(int, ConnectedVariable("squared.#state")) # should not require computation
        setstate!(int, ConnectedVariable("squared.#state"), 6.0)
        setstate!(int, ConnectedVariable("squared.#state"), 7.0)
        setstate!(int, ConnectedVariable("squared.#state"), 8.0)
        getstate(int, ConnectedVariable("squared.#state")) # required computation
        getstate(int) # should not require computation
        getstate(int, ConnectedVariable("squared.#state")) # should not require computation
        return nothing
    end

    # Make sure code is compiled before timing
    timer(int)
    int = init(comp)
    t = @timed timer(int)
    @test 0.2 < t.time < 0.3 # calling function f is required twice
end
