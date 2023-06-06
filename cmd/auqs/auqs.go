package main

import (
	"fmt"
	"math/cmplx"
	"math"
	"math/rand"
)

type Gates struct {
	i complex128
	singleQubitGates map[string][][]complex128
}

func NewGates() *Gates {
	g := &Gates{
		i: complex(0, 1),
		singleQubitGates: make(map[string][][]complex128),
	}

	// Pauli-X
	g.singleQubitGates["X"] = [][]complex128{
		{complex(0,0), complex(1,0)}, 
		{complex(1,0), complex(0,0)},
	}
	// Pauli-Y
	g.singleQubitGates["Y"] = [][]complex128{
		{complex(0,0), -g.i}, 
		{g.i, complex(0,0)},
	}
	// Pauli-Z
	g.singleQubitGates["Z"] = [][]complex128{
		{complex(1,0), complex(0,0)}, 
		{complex(0,0), complex(-1,0)},
	}
	// Hadamard
	h := 1 / math.Sqrt(2)
	g.singleQubitGates["H"] = [][]complex128{
		{complex(h, 0), complex(h,0)}, 
		{complex(h, 0), complex(-h,0)},
	}
	// Identity
	g.singleQubitGates["Id"] = [][]complex128{
		{complex(1,0), complex(0,0)}, 
		{complex(0,0), complex(1,0)},
	}
	// S
	g.singleQubitGates["S"] = [][]complex128{
		{complex(1, 0), complex(0, 0)},
		{complex(0, 0), g.i},
	}
	// Implement S Dagger
	// S Dagger
	g.singleQubitGates["SDagger"] = [][]complex128{
		{complex(1, 0), complex(0, 0)},
		{complex(0, 0), -g.i},
	}
	// T
	g.singleQubitGates["T"] = [][]complex128{
		{complex(1, 0), complex(0, 0)},
		{complex(0, 0), cmplx.Exp(g.i * math.Pi / 4)},
	}
	// TODO: implement T Dagger
	// T Dagger
	g.singleQubitGates["TDagger"] = [][]complex128{
		{complex(1, 0), complex(0, 0)},
		{complex(0, 0), cmplx.Exp(-g.i * math.Pi / 4,)},
	}

	return g
}

func GenerateGate(gate string, numQubits int, qubit1 int, qubit2 int) [][]complex128 {
	g := NewGates()

	if gate == "CNOT" {
		control := qubit1
		target := qubit2

		identity := g.singleQubitGates["Id"]
		X := g.singleQubitGates["X"]

		C := [][]complex128{
			{complex(math.NaN(), 0), complex(0, 0)}, 
			{complex(0, 0), complex(1, 0)},
		}

		gateOrder := make([][][]complex128, numQubits)
		for i := 0; i < numQubits; i++ {
			if i == control {
				gateOrder[i] = C
			} else if i == target {
				gateOrder[i] = X
			} else {
				gateOrder[i] = identity
			}
		}

		newGate := reduce(gateOrder)

		n := len(newGate)
		result := make([][]complex128, n)
		for i := 0; i < n; i++ {
			result[i] = make([]complex128, n)
			for j := 0; j < n; j++ {
				if !math.IsNaN(real(newGate[i][j])) && !math.IsNaN(imag(newGate[i][j])) {
					result[i][j] = newGate[i][j]
				} else if i == j {
					result[i][j] = complex(1, 0)
				} else {
					result[i][j] = complex(0, 0)
				}
			}
		}

		return result
	} else {
		identity := g.singleQubitGates["Id"]
		mainGate := g.singleQubitGates[gate]

		gateOrder := make([][][]complex128, numQubits)
		for i := 0; i < numQubits; i++ {
			if i == qubit1 {
				gateOrder[i] = mainGate
			} else {
				gateOrder[i] = identity
			}
		}

		return reduce(gateOrder)
	}
}

func reduce(gateOrder [][][]complex128) [][]complex128 {
	n := len(gateOrder)
	if n == 1 {
		return gateOrder[0]
	} 
	
	gate1 := gateOrder[0]
	gate2 := gateOrder[1]
	remaining := gateOrder[2:]
	
	result := kron(gate1, gate2)
	return reduce(append([][][]complex128{result}, remaining...))
}

func kron(a,b [][]complex128) [][]complex128 {
	n := len(a) * len(b)
	result := make([][]complex128, n)
	for i := 0; i < n; i++ {
		result[i] = make([]complex128, n)
		for j := 0; j < n; j++ {
			rowA := i / len(b)
			colA := j / len(b)
			rowB := i % len(b)
			colB := j % len(b)

			result[i][j] = a[rowA][colA] * b[rowB][colB]
		}
	}
	return result
}

type QuantumRegister struct {
	numQubits    int
	amplitudes   []complex128
	value        string
	probabilities []float64
}

func NewQuantumRegister(numQubits int) *QuantumRegister {
	qr := &QuantumRegister{
		numQubits:    numQubits,
		amplitudes:   make([]complex128, 1<<numQubits),
		value:        "",
		probabilities: nil,
	}

	// Set the probability of getting 0 when measured to 1
	qr.amplitudes[0] = complex(1, 0)

	return qr
}

func (qr *QuantumRegister) applyGate(gate string, qubit1 int, qubit2 int) {
	if qr.value != "" {
		panic("Cannot Apply Gate to Measured Register")
	} else {
		gateMatrix := GenerateGate(gate, qr.numQubits, qubit1, qubit2)
		qr.amplitudes = matrixMultiply(qr.amplitudes, gateMatrix)
	}
}

func (qr *QuantumRegister) measure() string {
	if qr.value != "" {
		return qr.value
	} else {
		qr.probabilities = make([]float64, len(qr.amplitudes))
		for i, amp := range qr.amplitudes {
			probability := math.Pow(cmplx.Abs(amp), 2)
			qr.probabilities[i] = probability
		}

		results := make([]int, len(qr.probabilities))
		for i := range results {
			results[i] = i
		}

		rand.Seed(42)

		index := rand.Intn(len(results))
		qr.value = fmt.Sprintf("%0*b", qr.numQubits, index)

		return qr.value
	}
}

func matrixMultiply(a []complex128, b [][]complex128) []complex128 {
	n := len(b)
	result := make([]complex128, n)
	for i := 0; i < n; i++ {
		result[i] = complex(0, 0)
		for j := 0; j < n; j++ {
			result[i] += a[j] * b[j][i]
		}
	}
	return result
}


func main() {
	fmt.Println("Execute Swap Algorithm")

	Swap := NewQuantumRegister(2)
	Swap.applyGate("X", 1, -1)
	
	// Value here is |10>
	// Start swapping
	Swap.applyGate("CNOT", 1, 2)
	Swap.applyGate("H", 1, -1)
	Swap.applyGate("H", 2, -1)
	Swap.applyGate("CNOT", 1, 2)
	Swap.applyGate("H", 1, -1)
	Swap.applyGate("H", 2, -1)
	Swap.applyGate("CNOT", 1, 2)  
	//End of swap
	fmt.Println("SWAP: |" + Swap.measure() + ">")
}