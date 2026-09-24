package optimizer

import (
	"math"
	"math/rand"
	"testing"
)

func TestGetExponentialTemperatureSteps(t *testing.T) {
	steps := GetExponentialTemperatureSteps(10.0, 5)
	if len(steps) != 5 {
		t.Fatalf("len(steps) = %d, want 5", len(steps))
	}
	if math.Abs(steps[0]-10.0) > 1e-9 {
		t.Errorf("steps[0] = %v, want 10.0 (the start temperature)", steps[0])
	}
	for i := 1; i < len(steps); i++ {
		if steps[i] >= steps[i-1] {
			t.Fatalf("steps[%d]=%v is not strictly less than steps[%d]=%v; schedule must decrease", i, steps[i], i-1, steps[i-1])
		}
	}
	if steps[len(steps)-1] >= steps[0] {
		t.Errorf("last step %v should be well below the start temperature %v", steps[len(steps)-1], steps[0])
	}
}

// A strictly better or equal move (log_p_next >= log_p_curr) must always be
// accepted, regardless of temperature or randomness.
func TestMoveAcceptanceAlwaysAcceptsBetterOrEqualMoves(t *testing.T) {
	random := rand.New(rand.NewSource(1))
	for _, temp := range []float64{0.0001, 1.0, 1000.0} {
		accept := getMoveAcceptanceFunc(temp, random)
		cases := []struct{ curr, next float64 }{
			{-5.0, -1.0}, // strictly better
			{-5.0, -5.0}, // equal
			{0.0, 0.0},   // equal
			{-100.0, 50.0},
		}
		for _, c := range cases {
			for range 10 { // randomness shouldn't matter here
				if !accept(c.curr, c.next) {
					t.Errorf("temp=%v: accept(%v, %v) = false, want true for a better-or-equal move", temp, c.curr, c.next)
				}
			}
		}
	}
}

// A worse move should be accepted far more often at high temperature than
// at low temperature (simulated annealing's defining property).
func TestMoveAcceptanceTemperatureEffect(t *testing.T) {
	const trials = 2000
	const curr, next = 0.0, -5.0 // next is strictly worse than curr

	acceptanceRate := func(temp float64, seed int64) float64 {
		random := rand.New(rand.NewSource(seed))
		accept := getMoveAcceptanceFunc(temp, random)
		accepted := 0
		for range trials {
			if accept(curr, next) {
				accepted++
			}
		}
		return float64(accepted) / trials
	}

	lowTempRate := acceptanceRate(0.01, 1)
	highTempRate := acceptanceRate(100.0, 2)

	if lowTempRate >= highTempRate {
		t.Errorf("low-temp acceptance rate %.3f should be well below high-temp rate %.3f", lowTempRate, highTempRate)
	}
	if lowTempRate > 0.05 {
		t.Errorf("low-temp acceptance rate %.3f too high; a much worse move should rarely be accepted when cold", lowTempRate)
	}
	if highTempRate < 0.5 {
		t.Errorf("high-temp acceptance rate %.3f too low; moves should be accepted often when hot", highTempRate)
	}
}

// fakeOptimizable counts how many times GenerateMove is invoked, and
// returns a new value each time so Optimize's state-threading is exercised.
type fakeOptimizable struct {
	calls *int
	value int
}

func (f fakeOptimizable) GenerateMove(accept_move func(log_p_curr, log_p_next float64) bool) Optimizable {
	*f.calls++
	// Always propose a strictly better move so it's deterministically accepted.
	accept_move(0.0, 1.0)
	return fakeOptimizable{calls: f.calls, value: f.value + 1}
}

func TestOptimizeCallsGenerateMoveExpectedNumberOfTimes(t *testing.T) {
	calls := 0
	start := fakeOptimizable{calls: &calls, value: 0}

	const temperatureSteps = 3
	const stepsPerTemp = 4

	result := Optimize(start, 10.0, temperatureSteps, stepsPerTemp)

	want := temperatureSteps * stepsPerTemp
	if calls != want {
		t.Errorf("GenerateMove called %d times, want %d (temperature_steps * steps_per_temp)", calls, want)
	}

	final, ok := result.(fakeOptimizable)
	if !ok {
		t.Fatalf("Optimize returned %T, want fakeOptimizable", result)
	}
	if final.value != want {
		t.Errorf("final value = %d, want %d; state was not threaded through every GenerateMove call", final.value, want)
	}
}
