#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// Dear ImGui
#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <string>
#include <algorithm>
#include <random>
#include <memory>
#include <cmath>

using namespace std;

// ------------------------------------------------------
// Configuration Constants
// ------------------------------------------------------
static const unsigned int WINDOW_WIDTH  = 960;
static const unsigned int WINDOW_HEIGHT = 600;
static const float BAR_GAP = 2.0f; // Gap between bars in pixels

// Default sorting parameters
static const int DEFAULT_NUMBER_OF_BARS = 100;
static chrono::milliseconds DEFAULT_SORT_DELAY(15);

// Names and details for each algorithm
static const char* ALGORITHM_NAMES[] = {
    "Bubble Sort",
    "Insertion Sort",
    "Quick Sort",
    "Merge Sort",
    "Heap Sort",
    "Selection Sort",
    "Shell Sort",
    "Cocktail Sort",
    "Comb Sort",
    "Odd-Even Sort"
};

static const char* ALGORITHM_DETAILS[] = {
    "Bubble Sort Detailed Steps:\n"
    "1. Compare adjacent elements and swap them if they are out of order.\n"
    "2. Continue comparing and swapping for the entire array.\n"
    "3. After each pass, the largest unsorted element bubbles up to its correct position.\n"
    "4. Repeat until no swaps are needed.",
    "Insertion Sort Detailed Steps:\n"
    "1. Start from the second element (the key).\n"
    "2. Compare the key with elements before it and shift larger elements to the right.\n"
    "3. Insert the key into its correct position.\n"
    "4. Repeat for all elements until the array is sorted.",
    "Quick Sort Detailed Steps:\n"
    "1. Choose a pivot (commonly the last element).\n"
    "2. Partition the array so that all elements less than the pivot come before it and all greater come after.\n"
    "3. Place the pivot in its final sorted position.\n"
    "4. Recursively apply the above steps to the sub-arrays on either side of the pivot.",
    "Merge Sort Detailed Steps:\n"
    "1. Divide the array into two roughly equal halves.\n"
    "2. Recursively sort each half.\n"
    "3. Merge the two sorted halves into one sorted array.\n"
    "4. Repeat until the entire array is merged into a sorted sequence.",
    "Heap Sort Detailed Steps:\n"
    "1. Build a max heap from the array.\n"
    "2. Swap the root (maximum element) with the last element of the heap.\n"
    "3. Reduce the heap size and heapify the new root to maintain the max heap property.\n"
    "4. Repeat until the heap size is 1 (sorted array).",
    "Selection Sort Detailed Steps:\n"
    "1. Find the minimum element in the unsorted portion of the array.\n"
    "2. Swap it with the first unsorted element.\n"
    "3. Move the boundary of the sorted portion one element to the right.\n"
    "4. Repeat until the entire array is sorted.",
    "Shell Sort Detailed Steps:\n"
    "1. Start with a gap (commonly half the array size).\n"
    "2. Perform a gapped insertion sort for the current gap value.\n"
    "3. Reduce the gap (often by half) and repeat the process.\n"
    "4. When the gap becomes 1, perform a final insertion sort.",
    "Cocktail Sort Detailed Steps:\n"
    "1. Traverse the array from left to right, swapping adjacent elements if necessary.\n"
    "2. Then traverse from right to left, again swapping where needed.\n"
    "3. Continue alternating forward and backward passes until no swaps occur.\n"
    "4. The array is sorted when a full cycle makes no swaps.",
    "Comb Sort Detailed Steps:\n"
    "1. Initialize the gap to the length of the array.\n"
    "2. Reduce the gap by a shrink factor (typically 1.3) until it becomes 1.\n"
    "3. Compare elements that are 'gap' apart and swap if they are in the wrong order.\n"
    "4. Repeat until a complete pass is made with no swaps.",
    "Odd-Even Sort Detailed Steps:\n"
    "1. In the even phase, compare and swap adjacent pairs starting from index 0.\n"
    "2. In the odd phase, do the same starting from index 1.\n"
    "3. Repeat both phases until no swaps occur in a complete cycle.\n"
    "4. The array is sorted when a full pass results in no swaps."
};

// ------------------------------------------------------
// Utility Function: Check if an array is sorted
// ------------------------------------------------------
static bool isArraySorted(const vector<int>& array) {
    for (size_t i = 0; i < array.size() - 1; i++) {
        if (array[i] > array[i + 1])
            return false;
    }
    return true;
}

// ------------------------------------------------------
// Base class for Sorting Algorithms (Abstract)
// ------------------------------------------------------
class Sorter {
public:
    // Constructor accepts the data to be sorted.
    Sorter(vector<int>& dataRef) : data(dataRef), comparisons(0), swaps(0), finished(false) {}
    virtual ~Sorter() {}

    // Perform one step of the sorting algorithm.
    virtual void step() = 0;

    // Should return true if sorting is complete.
    bool isFinished() const { return finished; }

    // Returns the number of comparisons and swaps.
    long long getComparisons() const { return comparisons; }
    long long getSwaps() const { return swaps; }

    // Returns the current state of the data.
    const vector<int>& getData() const { return data; }

    // Returns whether a given index is actively being operated on (to be highlighted).
    virtual bool isIndexHighlighted(int index) const = 0;

protected:
    vector<int>& data;  // Reference to the array being sorted.
    long long comparisons;
    long long swaps;
    bool finished;
};

// ------------------------------------------------------
// Bubble Sort Implementation
// ------------------------------------------------------
class BubbleSorter : public Sorter {
public:
    BubbleSorter(vector<int>& dataRef) 
        : Sorter(dataRef), i(0), j(0) {}

    void step() override {
        if (i >= static_cast<int>(data.size()) - 1) {
            finished = true;
            return;
        }
        comparisons++;
        if (data[j] > data[j + 1]) {
            swap(data[j], data[j + 1]);
            swaps++;
        }
        j++;
        if (j >= static_cast<int>(data.size()) - i - 1) {
            j = 0;
            i++;
        }
        if (isArraySorted(data))
            finished = true;
    }

    bool isIndexHighlighted(int index) const override {
        // Highlight the two adjacent indices being compared.
        return (index == j || index == j + 1);
    }

private:
    int i, j; // i: pass count, j: current index in pass.
};

// ------------------------------------------------------
// Insertion Sort Implementation
// ------------------------------------------------------
class InsertionSorter : public Sorter {
public:
    InsertionSorter(vector<int>& dataRef) 
        : Sorter(dataRef), i(1), j(0), key(0), keyHeld(false) {}

    void step() override {
        if (i >= static_cast<int>(data.size())) {
            finished = true;
            return;
        }
        if (!keyHeld) {
            key = data[i];
            j = i - 1;
            keyHeld = true;
        }
        if (j >= 0) comparisons++;
        if (j >= 0 && data[j] > key) {
            data[j + 1] = data[j];
            j--;
            swaps++;
        } else {
            data[j + 1] = key;
            i++;
            keyHeld = false;
            swaps++;
        }
    }

    bool isIndexHighlighted(int index) const override {
        // Highlight current key index and the index being compared.
        return (index == i || index == j);
    }

private:
    int i, j, key;
    bool keyHeld;
};

// ------------------------------------------------------
// Quick Sort Implementation
// ------------------------------------------------------
struct Range {
    int left, right;
};

class QuickSorter : public Sorter {
public:
    QuickSorter(vector<int>& dataRef) 
        : Sorter(dataRef), inPartition(false), leftBound(0), rightBound(0), pivotIndex(0)
    {
        stack.push_back({0, static_cast<int>(data.size()) - 1});
    }

    void step() override {
        if (!inPartition) {
            if (stack.empty()) {
                finished = true;
                return;
            }
            Range currentRange = stack.back();
            stack.pop_back();
            leftBound = currentRange.left;
            rightBound = currentRange.right;
            pivotIndex = rightBound;
            inPartition = true;
            return;
        }
        int p = partition(leftBound, rightBound);
        if (p - 1 > leftBound)
            stack.push_back({leftBound, p - 1});
        if (p + 1 < rightBound)
            stack.push_back({p + 1, rightBound});
        inPartition = false;
        if (isArraySorted(data))
            finished = true;
    }

    bool isIndexHighlighted(int index) const override {
        // Highlight the pivot element.
        return (index == pivotIndex);
    }

private:
    vector<Range> stack;
    bool inPartition;
    int leftBound, rightBound;
    int pivotIndex;

    int partition(int left, int right) {
        int pivotValue = data[right];
        int i = left - 1;
        for (int j = left; j < right; j++) {
            comparisons++;
            if (data[j] <= pivotValue) {
                i++;
                swap(data[i], data[j]);
                swaps++;
            }
        }
        swap(data[i+1], data[right]);
        swaps++;
        return i + 1;
    }
};

// ------------------------------------------------------
// Merge Sort Implementation
// ------------------------------------------------------
class MergeSorter : public Sorter {
public:
    MergeSorter(vector<int>& dataRef) 
        : Sorter(dataRef), subarraySize(1), inSubMerge(false), leftStart(0), mid(0), rightEnd(0)
    {
        mergeTemp.resize(data.size());
    }

    void step() override {
        if (subarraySize >= static_cast<int>(data.size())) {
            finished = true;
            return;
        }
        if (!inSubMerge) {
            if (leftStart >= static_cast<int>(data.size()) - 1) {
                subarraySize *= 2;
                leftStart = 0;
                if (subarraySize >= static_cast<int>(data.size()))
                    finished = true;
                return;
            }
            inSubMerge = true;
            mid = min(leftStart + subarraySize - 1, static_cast<int>(data.size()) - 1);
            rightEnd = min(leftStart + 2 * subarraySize - 1, static_cast<int>(data.size()) - 1);
            return;
        } else {
            mergeSubarray(leftStart, mid, rightEnd);
            inSubMerge = false;
            leftStart += 2 * subarraySize;
            if (isArraySorted(data))
                finished = true;
        }
    }

    bool isIndexHighlighted(int index) const override {
        // During merge, highlight indices in the current subarray being merged.
        if (inSubMerge && index >= leftStart && index <= rightEnd)
            return true;
        return false;
    }

private:
    int subarraySize;
    bool inSubMerge;
    int leftStart, mid, rightEnd;
    vector<int> mergeTemp;

    void mergeSubarray(int left, int mid, int right) {
        int i = left, j = mid + 1, k = left;
        while (i <= mid && j <= right) {
            comparisons++;
            if (data[i] <= data[j])
                mergeTemp[k++] = data[i++];
            else {
                mergeTemp[k++] = data[j++];
                swaps++;
            }
        }
        while (i <= mid)
            mergeTemp[k++] = data[i++];
        while (j <= right)
            mergeTemp[k++] = data[j++];
        for (int idx = left; idx <= right; idx++)
            data[idx] = mergeTemp[idx];
    }
};

// ------------------------------------------------------
// Heap Sort Implementation
// ------------------------------------------------------
class HeapSorter : public Sorter {
public:
    HeapSorter(vector<int>& dataRef)
        : Sorter(dataRef), heapBuilt(false), heapifyIndex(static_cast<int>(data.size()) / 2 - 1), heapSize(static_cast<int>(data.size()))
    {}

    void step() override {
        if (!heapBuilt) {
            if (heapifyIndex < 0) {
                heapBuilt = true;
                heapifyIndex = heapSize - 1;
                return;
            }
            heapify(heapSize, heapifyIndex);
            heapifyIndex--;
        } else {
            if (heapSize <= 1) {
                finished = true;
                return;
            }
            swap(data[0], data[heapSize - 1]);
            swaps++;
            heapSize--;
            heapify(heapSize, 0);
        }
    }

    bool isIndexHighlighted(int index) const override {
        // Always highlight the root (index 0) during heap sort.
        return (index == 0);
    }

private:
    bool heapBuilt;
    int heapifyIndex;
    int heapSize;

    void heapify(int n, int i) {
        int largest = i;
        int leftChild = 2 * i + 1;
        int rightChild = 2 * i + 2;
        if (leftChild < n) {
            comparisons++;
            if (data[leftChild] > data[largest])
                largest = leftChild;
        }
        if (rightChild < n) {
            comparisons++;
            if (data[rightChild] > data[largest])
                largest = rightChild;
        }
        if (largest != i) {
            swap(data[i], data[largest]);
            swaps++;
            heapify(n, largest);
        }
    }
};

// ------------------------------------------------------
// Selection Sort Implementation
// ------------------------------------------------------
class SelectionSorter : public Sorter {
public:
    SelectionSorter(vector<int>& dataRef)
        : Sorter(dataRef), i(0), j(0), minIndex(0)
    {}

    void step() override {
        if (i >= static_cast<int>(data.size()) - 1) {
            finished = true;
            return;
        }
        if (j < static_cast<int>(data.size())) {
            comparisons++;
            if (data[j] < data[minIndex])
                minIndex = j;
            j++;
            if (j >= static_cast<int>(data.size())) {
                if (minIndex != i) {
                    swap(data[i], data[minIndex]);
                    swaps++;
                }
                i++;
                j = i;
                minIndex = i;
            }
        }
    }

    bool isIndexHighlighted(int index) const override {
        // Highlight the current minimum index and the current comparison index.
        return (index == minIndex || index == j);
    }

private:
    int i, j, minIndex;
};

// ------------------------------------------------------
// Shell Sort Implementation
// ------------------------------------------------------
class ShellSorter : public Sorter {
public:
    ShellSorter(vector<int>& dataRef)
        : Sorter(dataRef), gap(static_cast<int>(data.size()) / 2), i(gap), inInnerLoop(false), currentJ(0), temp(0)
    {}

    void step() override {
        if (gap < 1) {
            finished = true;
            return;
        }
        if (!inInnerLoop) {
            if (i < static_cast<int>(data.size())) {
                temp = data[i];
                currentJ = i;
                inInnerLoop = true;
            } else {
                gap = gap / 2;
                if (gap < 1) {
                    finished = true;
                    return;
                }
                i = gap;
                return;
            }
        }
        if (inInnerLoop) {
            if (currentJ >= gap) {
                comparisons++;
                if (data[currentJ - gap] > temp) {
                    data[currentJ] = data[currentJ - gap];
                    swaps++;
                    currentJ -= gap;
                    return;
                }
            }
            data[currentJ] = temp;
            inInnerLoop = false;
            i++;
        }
    }

    bool isIndexHighlighted(int index) const override {
        // If not in inner loop, highlight current index; if in inner loop, highlight the active index and its gap partner.
        if (!inInnerLoop) {
            return (index == i);
        } else {
            return (index == currentJ || (currentJ - gap >= 0 && index == currentJ - gap));
        }
    }

private:
    int gap, i, currentJ, temp;
    bool inInnerLoop;
};

// ------------------------------------------------------
// Cocktail Sort Implementation
// ------------------------------------------------------
class CocktailSorter : public Sorter {
public:
    CocktailSorter(vector<int>& dataRef)
        : Sorter(dataRef), start(0), end(static_cast<int>(data.size()) - 1), currentIndex(0), forward(true)
    {
        currentIndex = start;
    }

    void step() override {
        if (start >= end) {
            finished = true;
            return;
        }
        if (forward) {
            if (currentIndex < end) {
                comparisons++;
                if (data[currentIndex] > data[currentIndex + 1]) {
                    swap(data[currentIndex], data[currentIndex + 1]);
                    swaps++;
                }
                currentIndex++;
            } else {
                forward = false;
                end--;
                currentIndex = end;
            }
        } else {
            if (currentIndex > start) {
                comparisons++;
                if (data[currentIndex - 1] > data[currentIndex]) {
                    swap(data[currentIndex - 1], data[currentIndex]);
                    swaps++;
                }
                currentIndex--;
            } else {
                forward = true;
                start++;
                currentIndex = start;
            }
        }
        if (isArraySorted(data))
            finished = true;
    }

    bool isIndexHighlighted(int index) const override {
        // Highlight the current index and its neighbor depending on the direction.
        if (forward)
            return (index == currentIndex || index == currentIndex + 1);
        else
            return (index == currentIndex || index == currentIndex - 1);
    }

private:
    int start, end, currentIndex;
    bool forward;
};

// ------------------------------------------------------
// Comb Sort Implementation
// ------------------------------------------------------
class CombSorter : public Sorter {
public:
    CombSorter(vector<int>& dataRef)
        : Sorter(dataRef), gap(static_cast<int>(data.size())), currentIndex(0), swapped(false)
    {}

    void step() override {
        if (currentIndex < static_cast<int>(data.size()) - gap) {
            comparisons++;
            if (data[currentIndex] > data[currentIndex + gap]) {
                swap(data[currentIndex], data[currentIndex + gap]);
                swaps++;
                swapped = true;
            }
            currentIndex++;
        } else {
            if (gap > 1) {
                gap = static_cast<int>(gap / 1.3f);
                if (gap < 1)
                    gap = 1;
            }
            if (gap == 1 && !swapped)
                finished = true;
            currentIndex = 0;
            swapped = false;
        }
    }

    bool isIndexHighlighted(int index) const override {
        // Highlight current index and the element gap away.
        return (index == currentIndex || index == currentIndex + gap);
    }

private:
    int gap, currentIndex;
    bool swapped;
};

// ------------------------------------------------------
// Odd-Even Sort Implementation
// ------------------------------------------------------
class OddEvenSorter : public Sorter {
public:
    OddEvenSorter(vector<int>& dataRef)
        : Sorter(dataRef), phase(0), indexOffset(0), swapped(false)
    {}

    void step() override {
        if (isArraySorted(data)) {
            finished = true;
            return;
        }
        int start = (phase == 0) ? 0 : 1;
        if (indexOffset < static_cast<int>(data.size()) - 1) {
            int idx = start + indexOffset;
            if (idx + 1 < static_cast<int>(data.size())) {
                comparisons++;
                if (data[idx] > data[idx + 1]) {
                    swap(data[idx], data[idx + 1]);
                    swaps++;
                    swapped = true;
                }
            }
            indexOffset++;
        } else {
            if (!swapped)
                finished = true;
            indexOffset = 0;
            phase = 1 - phase;
            swapped = false;
        }
    }

    bool isIndexHighlighted(int index) const override {
        int start = (phase == 0) ? 0 : 1;
        int currentIdx = start + indexOffset;
        return (index == currentIdx || index == currentIdx + 1);
    }

private:
    int phase;       // 0: even phase, 1: odd phase
    int indexOffset; // current offset in the phase
    bool swapped;
};

// ------------------------------------------------------
// SortingVisualizer: Manages the data, sorter and rendering
// ------------------------------------------------------
class SortingVisualizer {
public:
    SortingVisualizer()
        : window(nullptr),
          numberOfBars(DEFAULT_NUMBER_OF_BARS),
          sortDelay(DEFAULT_SORT_DELAY),
          sortingActive(false),
          paused(false),
          sortFinished(false),
          timeStopped(false),
          finalTime(0.0),
          currentAlgorithmIndex(0),
          finalPassIndex(0),
          lastFrameTime(0.0f),
          waveMode(false)
    {
        sortData.resize(numberOfBars);
        displayHeights.resize(numberOfBars, 0.0f);
    }

    ~SortingVisualizer() {
        // Cleanup ImGui and GLFW
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glDeleteProgram(shaderProgram);
        glDeleteVertexArrays(1, &barVAO);
        glDeleteBuffers(1, &barVBO);
        glfwTerminate();
    }

    // Initialize GLFW, GLAD, ImGui and OpenGL objects.
    bool initialize() {
        if (!glfwInit()) {
            cerr << "Failed to initialize GLFW\n";
            return false;
        }
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Sorting Visualization", nullptr, nullptr);
        if (!window) {
            cerr << "Failed to create GLFW window\n";
            glfwTerminate();
            return false;
        }
        glfwMakeContextCurrent(window);
        glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);

        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
            cerr << "Failed to initialize GLAD\n";
            return false;
        }

        // Initialize ImGui
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        setCustomTheme();
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 330 core");

        // Setup OpenGL objects (shaders, VAO, VBO)
        setupShaders();
        setupBarQuad();

        // Initialize default data and sorter
        randomizeData();
        resetSorterState();

        // Initialize the frame time tracker for smooth interpolation.
        lastFrameTime = static_cast<float>(glfwGetTime());

        return true;
    }

    // Main loop: handles events, sorting steps and rendering.
    void run() {
        bool spacePreviouslyPressed = false;
        sortStartTime = chrono::steady_clock::now();
        lastStepTime = sortStartTime;
        finalPassLastUpdate = sortStartTime;

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();

            // Check for space key to toggle pause/resume.
            bool spacePressed = (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS);
            if (spacePressed && !spacePreviouslyPressed && sortingActive)
                paused = !paused;
            spacePreviouslyPressed = spacePressed;

            // Execute sorting step if active, not paused and not finished.
            if (sortingActive && !paused && !sortFinished && !waveMode) {
                auto now = chrono::steady_clock::now();
                if (chrono::duration_cast<chrono::milliseconds>(now - lastStepTime) >= sortDelay) {
                    if (currentSorter)
                        currentSorter->step();
                    lastStepTime = now;
                    if (currentSorter && currentSorter->isFinished() && !timeStopped) {
                        finalTime = chrono::duration_cast<chrono::milliseconds>(now - sortStartTime).count() / 1000.0;
                        timeStopped = true;
                        sortFinished = true;
                    }
                }
            }

            // Update display heights for smooth animation or wave mode.
            updateDisplayHeights();

            // Handle final pass animation after sort is complete.
            if (sortFinished && finalPassIndex < numberOfBars) {
                auto now = chrono::steady_clock::now();
                // Increase finalPassIndex periodically
                if (chrono::duration_cast<chrono::milliseconds>(now - finalPassLastUpdate).count() > 20) {
                    finalPassIndex++;
                    finalPassLastUpdate = now;
                }
            }

            // Render the ImGui interface and bars.
            render();
            glfwSwapBuffers(window);
        }
    }

private:
    GLFWwindow* window;
    int numberOfBars;
    vector<int> sortData;         // Data to be sorted.
    vector<float> displayHeights; // Current heights used for bar animation.
    unique_ptr<Sorter> currentSorter; // Polymorphic sorter instance.
    chrono::steady_clock::time_point sortStartTime, lastStepTime, finalPassLastUpdate;
    chrono::milliseconds sortDelay;
    bool sortingActive, paused, sortFinished, timeStopped;
    double finalTime;
    int currentAlgorithmIndex; // Selected algorithm index (from ImGui combo).
    int finalPassIndex;        // Used for final pass animation.
    float lastFrameTime;       // Tracks the last frame time for smooth interpolation.

    // --- Easter Egg State ---
    bool waveMode;                  // When true, show wave animation.
    vector<int> buttonCombo;        // Records the order of button presses.

    // OpenGL objects for rendering bars.
    unsigned int shaderProgram, barVAO, barVBO;

    // --------------------------------------------------
    // Data Initialization Methods
    // --------------------------------------------------

    // Randomize data: fill with sorted numbers and then shuffle.
    void randomizeData() {
        sortData.resize(numberOfBars);
        displayHeights.resize(numberOfBars, 0.0f);
        for (int i = 0; i < numberOfBars; i++) {
            sortData[i] = i + 1;
        }
        random_device rd;
        mt19937 gen(rd());
        shuffle(sortData.begin(), sortData.end(), gen);
    }

    // Create a mountain pattern.
    void mountainData() {
        sortData.resize(numberOfBars);
        displayHeights.resize(numberOfBars, 0.0f);
        int mid = numberOfBars / 2;
        for (int i = 0; i < numberOfBars; i++) {
            if (i <= mid)
                sortData[i] = i + 1;
            else
                sortData[i] = sortData[mid - (i - mid)];
        }
        // Rescale to span from 1 to numberOfBars.
        int currentMin = *min_element(sortData.begin(), sortData.end());
        int currentMax = *max_element(sortData.begin(), sortData.end());
        for (int i = 0; i < numberOfBars; i++) {
            sortData[i] = 1 + (sortData[i] - currentMin) * (numberOfBars - 1) / (currentMax - currentMin);
        }
    }

    // Create a valley pattern.
    void valleyData() {
        sortData.resize(numberOfBars);
        displayHeights.resize(numberOfBars, 0.0f);
        int mid = numberOfBars / 2;
        for (int i = 0; i < numberOfBars; i++) {
            if (i <= mid)
                sortData[i] = mid - i + 1;
            else
                sortData[i] = i - mid + 1;
        }
        int currentMin = *min_element(sortData.begin(), sortData.end());
        int currentMax = *max_element(sortData.begin(), sortData.end());
        for (int i = 0; i < numberOfBars; i++) {
            sortData[i] = 1 + (sortData[i] - currentMin) * (numberOfBars - 1) / (currentMax - currentMin);
        }
    }

    // Create reverse-sorted data.
    void reverseData() {
        sortData.resize(numberOfBars);
        displayHeights.resize(numberOfBars, 0.0f);
        for (int i = 0; i < numberOfBars; i++) {
            sortData[i] = numberOfBars - i;
        }
    }

    // --------------------------------------------------
    // Sorter and State Reset
    // --------------------------------------------------
    void resetSorterState() {
        sortFinished = false;
        timeStopped = false;
        finalTime = 0.0;
        paused = false;
        finalPassIndex = 0;
        sortStartTime = chrono::steady_clock::now();
        lastStepTime = sortStartTime;
        finalPassLastUpdate = sortStartTime;
        // Create a new sorter based on the selected algorithm.
        createSorter(currentAlgorithmIndex);
    }

    void createSorter(int algoIndex) {
        // Use polymorphism to instantiate the correct sorter.
        switch(algoIndex) {
            case 0: currentSorter = make_unique<BubbleSorter>(sortData); break;
            case 1: currentSorter = make_unique<InsertionSorter>(sortData); break;
            case 2: currentSorter = make_unique<QuickSorter>(sortData); break;
            case 3: currentSorter = make_unique<MergeSorter>(sortData); break;
            case 4: currentSorter = make_unique<HeapSorter>(sortData); break;
            case 5: currentSorter = make_unique<SelectionSorter>(sortData); break;
            case 6: currentSorter = make_unique<ShellSorter>(sortData); break;
            case 7: currentSorter = make_unique<CocktailSorter>(sortData); break;
            case 8: currentSorter = make_unique<CombSorter>(sortData); break;
            case 9: currentSorter = make_unique<OddEvenSorter>(sortData); break;
            default: currentSorter = make_unique<BubbleSorter>(sortData); break;
        }
    }

    // --------------------------------------------------
    // Easter Egg: Record Button Presses and Check Combo
    // --------------------------------------------------
    void recordButtonPress(int code) {
        buttonCombo.push_back(code);
        // Keep only the last 6 button presses
        if(buttonCombo.size() > 6)
            buttonCombo.erase(buttonCombo.begin());
        // Check for combo: Randomize(0), Mountain(1), Mountain(1), Valley(2), Reverse(3), Valley(2)
        if(buttonCombo.size() == 6 &&
           buttonCombo[0] == 0 &&
           buttonCombo[1] == 1 &&
           buttonCombo[2] == 1 &&
           buttonCombo[3] == 2 &&
           buttonCombo[4] == 3 &&
           buttonCombo[5] == 2) {
            waveMode = true;
            // Clear the combo after triggering wave mode.
            buttonCombo.clear();
        }
    }

    // --------------------------------------------------
    // OpenGL and ImGui Rendering Helpers
    // --------------------------------------------------

    // GLFW framebuffer resize callback.
    static void framebufferSizeCallback(GLFWwindow* window, int width, int height) {
        glViewport(0, 0, width, height);
    }

    // Compile a shader from source.
    unsigned int compileShader(unsigned int type, const char* source) {
        unsigned int shader = glCreateShader(type);
        glShaderSource(shader, 1, &source, nullptr);
        glCompileShader(shader);
        int success;
        char infoLog[512];
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if(!success) {
            glGetShaderInfoLog(shader, 512, nullptr, infoLog);
            cerr << "Shader compile error:\n" << infoLog << endl;
        }
        return shader;
    }

    // Setup the shader program for drawing bars.
    void setupShaders() {
        const char* vertexShaderSource = R"(
        #version 330 core
        layout(location=0) in vec2 aPos;
        uniform mat4 uMVP;
        void main(){
            gl_Position = uMVP * vec4(aPos, 0.0, 1.0);
        }
        )";
        const char* fragmentShaderSource = R"(
        #version 330 core
        out vec4 FragColor;
        uniform vec4 uColor;
        void main(){
            FragColor = uColor;
        }
        )";
        unsigned int vertexShader = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
        unsigned int fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);
        shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);
        int success;
        char infoLog[512];
        glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
        if(!success) {
            glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog);
            cerr << "Shader program link error:\n" << infoLog << endl;
        }
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
    }

    // Setup the VAO and VBO for a single bar quad.
    void setupBarQuad() {
        float barQuadVertices[] = {
            0.0f, 0.0f,
            1.0f, 0.0f,
            1.0f, 1.0f,
            0.0f, 0.0f,
            1.0f, 1.0f,
            0.0f, 1.0f
        };
        glGenVertexArrays(1, &barVAO);
        glGenBuffers(1, &barVBO);
        glBindVertexArray(barVAO);
        glBindBuffer(GL_ARRAY_BUFFER, barVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(barQuadVertices), barQuadVertices, GL_STATIC_DRAW);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glBindVertexArray(0);
    }

    // Set a custom dark theme for ImGui.
    void setCustomTheme() {
        ImGuiStyle& style = ImGui::GetStyle();
        ImVec4* colors = style.Colors;
        colors[ImGuiCol_WindowBg]         = ImVec4(0.10f, 0.12f, 0.15f, 1.00f);
        colors[ImGuiCol_ChildBg]          = ImVec4(0.15f, 0.17f, 0.20f, 1.00f);
        colors[ImGuiCol_FrameBg]          = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
        colors[ImGuiCol_FrameBgHovered]   = ImVec4(0.30f, 0.32f, 0.37f, 1.00f);
        colors[ImGuiCol_FrameBgActive]    = ImVec4(0.25f, 0.27f, 0.32f, 1.00f);
        colors[ImGuiCol_TitleBg]          = ImVec4(0.12f, 0.15f, 0.20f, 1.00f);
        colors[ImGuiCol_TitleBgActive]    = ImVec4(0.15f, 0.18f, 0.23f, 1.00f);
        colors[ImGuiCol_Button]           = ImVec4(0.20f, 0.50f, 0.80f, 1.00f);
        colors[ImGuiCol_ButtonHovered]    = ImVec4(0.25f, 0.60f, 0.90f, 1.00f);
        colors[ImGuiCol_ButtonActive]     = ImVec4(0.15f, 0.45f, 0.75f, 1.00f);
        colors[ImGuiCol_Header]           = ImVec4(0.20f, 0.50f, 0.80f, 1.00f);
        colors[ImGuiCol_HeaderHovered]    = ImVec4(0.25f, 0.60f, 0.90f, 1.00f);
        colors[ImGuiCol_HeaderActive]     = ImVec4(0.15f, 0.45f, 0.75f, 1.00f);
        colors[ImGuiCol_SliderGrab]       = ImVec4(0.20f, 0.50f, 0.80f, 1.00f);
        colors[ImGuiCol_SliderGrabActive] = ImVec4(0.25f, 0.60f, 0.90f, 1.00f);
        colors[ImGuiCol_Tab]              = ImVec4(0.12f, 0.15f, 0.20f, 1.00f);
        colors[ImGuiCol_TabHovered]       = ImVec4(0.25f, 0.60f, 0.90f, 1.00f);
        colors[ImGuiCol_TabActive]        = ImVec4(0.20f, 0.50f, 0.80f, 1.00f);
    }

    // --------------------------------------------------
    // Update display heights based on current sort data or wave mode.
    // --------------------------------------------------
    void updateDisplayHeights() {
        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        float minHeight = 10.0f;
        float maxHeight = fbHeight * 0.9f;
        float currentTime = static_cast<float>(glfwGetTime());

        // If wave mode is active, use a sine wave for the bar heights.
        if(waveMode) {
            float base = (minHeight + maxHeight) / 2.0f;
            float amplitude = (maxHeight - minHeight) / 4.0f;
            float speed = 2.0f;
            float phaseOffset = 0.5f;
            for (int i = 0; i < numberOfBars; i++) {
                displayHeights[i] = base + amplitude * sin(currentTime * speed + i * phaseOffset);
            }
            return;
        }

        // Otherwise, smoothly animate from the current displayHeights to the target.
        float dt = currentTime - lastFrameTime;
        lastFrameTime = currentTime;
        float interpolationFactor = 1 - exp(-10.0f * dt);

        for (int i = 0; i < numberOfBars; i++) {
            float targetHeight = minHeight + (maxHeight - minHeight) * ((sortData[i] - 1) / float(numberOfBars - 1));
            displayHeights[i] = glm::mix(displayHeights[i], targetHeight, interpolationFactor);
        }
    }

    // --------------------------------------------------
    // Render the GUI and the sorting bars.
    // --------------------------------------------------
    void render() {
        // Start new ImGui frame.
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // ImGui Window with Tabs for Controls, Settings, Statistics, and About.
        ImGui::Begin("Sorting Visualizer");
        if (ImGui::BeginTabBar("MainTabBar")) {
            if (ImGui::BeginTabItem("Controls")) {
                // If wave mode is active, show a special UI.
                if(waveMode) {
                    ImGui::Text("Wave Mode Activated!");
                    if(ImGui::Button("Exit Wave Mode")) {
                        waveMode = false;
                        buttonCombo.clear();
                        randomizeData();
                        resetSorterState();
                        sortingActive = false;
                    }
                } else if (!sortingActive) {
                    // Select sorting algorithm.
                    if (ImGui::BeginCombo("Algorithm", ALGORITHM_NAMES[currentAlgorithmIndex])) {
                        for (int n = 0; n < IM_ARRAYSIZE(ALGORITHM_NAMES); n++) {
                            bool isSelected = (currentAlgorithmIndex == n);
                            if (ImGui::Selectable(ALGORITHM_NAMES[n], isSelected))
                                currentAlgorithmIndex = n;
                            if (isSelected)
                                ImGui::SetItemDefaultFocus();
                        }
                        ImGui::EndCombo();
                    }
                    // Pattern buttons
                    if (ImGui::Button("Randomize")) {
                        randomizeData();
                        resetSorterState();
                        recordButtonPress(0); // Code 0 for Randomize
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Mountain")) {
                        mountainData();
                        resetSorterState();
                        recordButtonPress(1); // Code 1 for Mountain
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Valley")) {
                        valleyData();
                        resetSorterState();
                        recordButtonPress(2); // Code 2 for Valley
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Reverse")) {
                        reverseData();
                        resetSorterState();
                        recordButtonPress(3); // Code 3 for Reverse
                    }
                    // Start sorting
                    if (ImGui::Button("Start")) {
                        createSorter(currentAlgorithmIndex);
                        resetSorterState();
                        sortingActive = true;
                        paused = false;
                        sortStartTime = chrono::steady_clock::now();
                        lastStepTime = sortStartTime;
                    }
                } else {
                    // Sorting is active
                    if (!paused) {
                        if (ImGui::Button("Pause"))
                            paused = true;
                    } else {
                        if (ImGui::Button("Resume")) {
                            paused = false;
                            lastStepTime = chrono::steady_clock::now();
                        }
                        ImGui::SameLine();
                        if (ImGui::Button("Step"))
                            if (currentSorter) currentSorter->step();
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Reset")) {
                        resetSorterState();
                        sortingActive = false;
                        paused = false;
                    }
                }
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Choose a pattern and start the sort");
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Settings")) {
                static float delayMs = float(sortDelay.count());
                if (ImGui::SliderFloat("Sort Speed (ms)", &delayMs, 1.0f, 100.0f, "%.0f"))
                    sortDelay = chrono::milliseconds(int(delayMs));
                static int bars = numberOfBars;
                if (ImGui::SliderInt("Number of Bars", &bars, 10, 300)) {
                    if (bars != numberOfBars) {
                        numberOfBars = bars;
                        randomizeData();
                        resetSorterState();
                    }
                }
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Statistics")) {
                ImGui::Text("Algorithm: %s", ALGORITHM_NAMES[currentAlgorithmIndex]);
                double elapsedSeconds = 0.0;
                if (sortFinished && timeStopped)
                    elapsedSeconds = finalTime;
                else if (sortingActive && !paused) {
                    auto now = chrono::steady_clock::now();
                    elapsedSeconds = chrono::duration_cast<chrono::milliseconds>(now - sortStartTime).count() / 1000.0;
                }
                ImGui::Text("Time: %.2f s", elapsedSeconds);
                if (currentSorter) {
                    ImGui::Text("Comparisons: %lld", currentSorter->getComparisons());
                    ImGui::Text("Swaps:       %lld", currentSorter->getSwaps());
                }
                float progress = 0.0f;
                ImGui::ProgressBar(progress, ImVec2(-1.0f, 0.0f));
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("About")) {
                ImGui::BeginChild("AboutChild", ImVec2(400, 300), true, ImGuiWindowFlags_HorizontalScrollbar);
                ImGui::TextWrapped("%s", ALGORITHM_DETAILS[currentAlgorithmIndex]);
                ImGui::Spacing();
                ImGui::TextWrapped("Time Complexity:\n"
                                   "Bubble Sort: O(n^2)\n"
                                   "Insertion Sort: O(n^2) (O(n) best-case)\n"
                                   "Quick Sort: O(n log n) average\n"
                                   "Merge Sort: O(n log n)\n"
                                   "Heap Sort: O(n log n)\n"
                                   "Selection Sort: O(n^2)\n"
                                   "Shell Sort: ~O(n log^2 n)\n"
                                   "Cocktail Sort: O(n^2)\n"
                                   "Comb Sort: O(n^2) worst-case\n"
                                   "Odd-Even Sort: O(n^2)");
                ImGui::EndChild();
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
        ImGui::End();

        // Clear screen.
        glClearColor(0.15f, 0.15f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Setup projection for 2D rendering.
        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        glm::mat4 projection = glm::ortho(0.0f, float(fbWidth), 0.0f, float(fbHeight), -1.0f, 1.0f);
        glUseProgram(shaderProgram);
        glBindVertexArray(barVAO);
        int locMVP = glGetUniformLocation(shaderProgram, "uMVP");
        int locColor = glGetUniformLocation(shaderProgram, "uColor");

        // Get mouse position for hover effect.
        double mouseX, mouseY;
        glfwGetCursorPos(window, &mouseX, &mouseY);
        float oglMouseY = float(fbHeight) - float(mouseY);
        float barWidth = (fbWidth - (numberOfBars + 1) * BAR_GAP) / float(numberOfBars);

        // Draw each bar.
        for (int i = 0; i < numberOfBars; i++) {
            float xPos = BAR_GAP + i * (barWidth + BAR_GAP);
            float barHeight = displayHeights[i];

            // Draw shadow.
            {
                glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(xPos + 2.0f, 2.0f, 0.0f));
                model = glm::scale(model, glm::vec3(barWidth, barHeight, 1.0f));
                glm::mat4 mvp = projection * model;
                glUniformMatrix4fv(locMVP, 1, GL_FALSE, glm::value_ptr(mvp));
                glUniform4f(locColor, 0.1f, 0.1f, 0.1f, 0.5f);
                glDrawArrays(GL_TRIANGLES, 0, 6);
            }

            // Determine bar color.
            glm::vec4 barColor(0.9f, 0.9f, 0.9f, 1.0f);

            // If sorting not finished, highlight current indices in red.
            if (!sortFinished) {
                if (currentSorter && currentSorter->isIndexHighlighted(i)) {
                    barColor = glm::vec4(0.8f, 0.2f, 0.2f, 1.0f);
                }
            }
            else {
                // Sorting finished: bars are white unless in final pass region.
                barColor = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);

                // If i <= finalPassIndex, apply a pulsing lighter green.
                if (i <= finalPassIndex) {
                    // Create a pulsing effect
                    float t = float(glfwGetTime()) * 5.0f;  // Adjust speed as needed
                    float pulse = 0.5f + 0.5f * sin(t);

                    // Base lighter green color
                    glm::vec4 baseColor(0.4f, 0.9f, 0.4f, 1.0f);

                    // We'll vary brightness slightly between 80% and 100%
                    float brightnessFactor = 0.8f + 0.2f * pulse;

                    barColor = glm::vec4(baseColor.r * brightnessFactor,
                                         baseColor.g * brightnessFactor,
                                         baseColor.b * brightnessFactor,
                                         1.0f);
                }
            }

            // Apply hover effect (slight yellow tint) on top of current color.
            bool isHover = (mouseX >= xPos && mouseX <= xPos + barWidth && oglMouseY <= barHeight);
            if (isHover) {
                barColor = glm::mix(barColor, glm::vec4(1.0f, 1.0f, 0.7f, 1.0f), 0.5f);
            }

            glUniform4fv(locColor, 1, &barColor[0]);

            // Set model transform and draw the bar.
            glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(xPos, 0.0f, 0.0f));
            model = glm::scale(model, glm::vec3(barWidth, barHeight, 1.0f));
            glm::mat4 mvp = projection * model;
            glUniformMatrix4fv(locMVP, 1, GL_FALSE, glm::value_ptr(mvp));
            glDrawArrays(GL_TRIANGLES, 0, 6);
        }

        // Render ImGui on top.
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }
};

// ------------------------------------------------------
// Main Entry Point
// ------------------------------------------------------
int main() {
    SortingVisualizer visualizer;
    if (!visualizer.initialize())
        return -1;
    visualizer.run();
    return 0;
}