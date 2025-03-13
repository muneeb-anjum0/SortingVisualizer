# Sorting Algorithm Visualizer

## Overview
Sorting Algorithm Visualizer is an interactive application that demonstrates various sorting algorithms in real-time using **OpenGL**, **GLFW**, **GLM**, and **Dear ImGui**. The project enables users to select different sorting algorithms, adjust sorting speed, and analyze their performance.

This tool serves as an **educational aid** for students, programmers, and algorithm enthusiasts to understand sorting techniques through **dynamic visual representation** and **Object-Oriented Programming (OOP)** concepts.
![Demo](https://github.com/user-attachments/assets/a403f7a5-5670-4373-ba2c-816ba87711b4)

## Features
- **Real-time visualization** of sorting algorithms.
- **Object-Oriented Design** for modular and maintainable code.
- **Interactive GUI** using **Dear ImGui**.
- **10 Popular Sorting Algorithms**:
  - Bubble Sort
  - Insertion Sort
  - Quick Sort
  - Merge Sort
  - Heap Sort
  - Selection Sort
  - Shell Sort
  - Cocktail Sort
  - Comb Sort
  - Odd-Even Sort
- **Encapsulation of Sorting Algorithms** into reusable classes.
- **Polymorphism** to allow dynamic selection of sorting algorithms.
- **Smooth animations** with interpolation.
- **Real-time statistics** (comparisons, swaps, elapsed time).
- **Customizable Parameters**:
  - Sorting speed
  - Number of elements
  - Initial array patterns (Random, Reverse, Mountain, Valley)
- **Hover to highlight bars** for better tracking.
- **Easter Egg Wave Mode** (Secret feature).

## Technologies Used
- **C++ (Object-Oriented Programming)**
- **OpenGL** (Rendering)
- **GLFW** (Window & input management)
- **GLAD** (OpenGL loader)
- **GLM** (Matrix & vector math)
- **Dear ImGui** (GUI & Controls)
- **STL Containers** (Vectors, Algorithms, Chrono, Random)


## Object-Oriented Programming (OOP) Concepts Used
The project follows **OOP principles** to make the code maintainable and scalable:

1. **Encapsulation**  
   - Sorting algorithms are implemented as separate classes.
   - The `Sorter` base class encapsulates common sorting behavior.
   - Each sorting algorithm (Bubble Sort, Quick Sort, etc.) is encapsulated in its own derived class.

2. **Inheritance**  
   - The `Sorter` class acts as a **base class**.
   - **Derived classes** (e.g., `BubbleSorter`, `QuickSorter`, `MergeSorter`) inherit from `Sorter`.

3. **Polymorphism**  
   - The project uses **virtual functions** to allow dynamic binding.
   - The `step()` method is overridden in each derived sorting class.
   - This allows the visualizer to **switch sorting algorithms dynamically**.

4. **Abstraction**  
   - The `Sorter` class defines an **abstract interface** for sorting.
   - Sorting implementations are hidden behind this interface.


## Installation & Running the Project
### Prerequisites
Ensure you have the following installed:
- **g++ (C++17 or later)**
- **OpenGL 3.3+**
- **GLFW** (`glfw3.h`)
- **GLAD** (`glad.h`)
- **GLM** (`glm.hpp`)
- **Dear ImGui** (`imgui.h`)

### Manually Compiling & Running (Without CMake)

#### **1. Install Dependencies**
Before compiling, ensure you have installed:
- **Linux (Debian-based)** (Ubuntu, Debian, etc.):
  sudo apt update
  sudo apt install g++ cmake xorg-dev libglfw3 libglfw3-dev libglm-dev
  
- **Windows**:  
  - Install [MinGW-w64](https://www.mingw-w64.org/)
  - Download [GLFW](https://www.glfw.org/)
  - Download [GLAD](https://glad.dav1d.de/)
  - Include these libraries in your project directory.

#### **2. Compile the Project**
Use the following `g++` command to compile:
g++ -std=c++17 -Wall -Wextra -O2 main.cpp -o SortingVisualizer -lglfw -lGL -ldl -lX11 -pthread

- `-std=c++17` → Uses C++17 standard.
- `-Wall -Wextra` → Enables warnings.
- `-O2` → Optimized compilation.
- `-lglfw` → Links **GLFW**.
- `-lGL -ldl -lX11 -pthread` → Required for **OpenGL**.

#### **3. Run the Executable**
./SortingVisualizer

## How to Use
1. **Choose a Sorting Algorithm** from the dropdown.
2. **Modify settings**:
   - Adjust **sorting speed**.
   - Change the **number of elements**.
   - Choose an **initial pattern** (Random, Reverse, Mountain, Valley).
3. **Start the sorting process**:
   - **Start**: Begins sorting animation.
   - **Pause/Resume**: Toggle with **Spacebar**.
   - **Step**: Manually progress sorting one step at a time.
   - **Reset**: Reshuffle and reset sorting.
4. **Analyze real-time statistics** (comparisons, swaps, elapsed time).

---

## Sorting Algorithms & Time Complexity

| Algorithm       | Best Case  | Average Case | Worst Case  |
|-----------------|------------|--------------|-------------|
| Bubble Sort     | O(n)       | O(n²)        | O(n²)       |
| Insertion Sort  | O(n)       | O(n²)        | O(n²)       |
| Quick Sort      | O(n log n) | O(n log n)   | O(n²) (rare)|
| Merge Sort      | O(n log n) | O(n log n)   | O(n log n)  |
| Heap Sort       | O(n log n) | O(n log n)   | O(n log n)  |
| Selection Sort  | O(n²)      | O(n²)        | O(n²)       |
| Shell Sort      | O(n log n) | O(n log² n)  | O(n²)       |
| Cocktail Sort   | O(n)       | O(n²)        | O(n²)       |
| Comb Sort       | O(n log n) | O(n log n)   | O(n²)       |
| Odd-Even Sort   | O(n)       | O(n²)        | O(n²)       |

## Future Improvements
- Add **Radix Sort & Counting Sort** for non-comparative sorting.
- Add **Parallel Sorting using Multithreading**.
- Add **Save & Load Sorting Animations as GIFs**.
- Improve **visual customization (themes, colors, UI tweaks)**.

## Contributing
Contributions are welcome! If you’d like to improve this project:
1. **Fork the repository**.
2. Create a **feature branch** (`git checkout -b feature-name`).
3. **Commit your changes** (`git commit -m "Added new feature"`).
4. **Push to the branch** (`git push origin feature-name`).
5. Open a **Pull Request**.

---

## License
This project is licensed under the **MIT License**.

---

## Acknowledgments
- **OpenGL & GLFW** for providing the foundation for rendering.
- **Dear ImGui** for the GUI interface.
- **Computer Science Community** for sorting algorithm research.
- Inspired by various **sorting visualizations** across the web.

---

## Contact
If you have any questions or suggestions, feel free to reach out:
- **Email:** muneeb.anjum@hotmail.com
- **GitHub:** https://github.com/muneeb-anjum0
- **LinkedIn:** https://www.linkedin.com/in/muneeb-anjum-a9386625a/

---
