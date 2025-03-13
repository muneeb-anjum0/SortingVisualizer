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

// Window and layout settings
static const unsigned int kWindowWidth  = 960;
static const unsigned int kWindowHeight = 600;
static const float kBarGap = 2.0f;

// Default sort settings
static const int kDefaultBarCount = 100;
static chrono::milliseconds kDefaultSortDelay(15);

// Names and info for each sorting algorithm
static const char* algorithmNames[] = {
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

static const char* algorithmDescriptions[] = {
    "Bubble Sort:\n1. Compare neighboring numbers and swap if needed.\n2. Repeat for the entire list until sorted.",
    "Insertion Sort:\n1. Start from the second number.\n2. Shift larger numbers right and insert the key in place.",
    "Quick Sort:\n1. Pick a pivot.\n2. Rearrange numbers around the pivot.\n3. Recursively sort the partitions.",
    "Merge Sort:\n1. Split the list into halves.\n2. Sort each half and merge them back together.",
    "Heap Sort:\n1. Build a max heap from the numbers.\n2. Swap the root with the last element, reduce heap size, and heapify.",
    "Selection Sort:\n1. Find the smallest number in the unsorted part.\n2. Swap it with the first unsorted number.",
    "Shell Sort:\n1. Use a gap to compare elements far apart.\n2. Reduce the gap until it becomes 1 and finish with insertion sort.",
    "Cocktail Sort:\n1. Traverse forward swapping when needed, then backward.\n2. Continue until no swaps are required.",
    "Comb Sort:\n1. Start with a large gap and reduce it gradually.\n2. Compare and swap distant numbers until the list is sorted.",
    "Odd-Even Sort:\n1. Alternate between comparing odd-even and even-odd pairs.\n2. Keep doing this until the list is in order."
};

// A simple helper to check if a vector is already sorted.
static bool isSorted(const vector<int>& nums) {
    for (size_t i = 0; i < nums.size() - 1; i++) {
        if (nums[i] > nums[i + 1])
            return false;
    }
    return true;
}

// ---------------------------------------------------------------------
// Abstract base class for our sorting algorithms.
// ---------------------------------------------------------------------
class SortAlgorithm {
public:
    SortAlgorithm(vector<int>& nums) 
        : numbers(nums), numComparisons(0), numSwaps(0), done(false) {}
    virtual ~SortAlgorithm() {}

    virtual void step() = 0;

    bool isDone() const { return done; }
    long long getComparisons() const { return numComparisons; }
    long long getSwaps() const { return numSwaps; }
    const vector<int>& getNumbers() const { return numbers; }

    // This helps highlight which indices are being worked on.
    virtual bool isIndexActive(int idx) const = 0;

protected:
    vector<int>& numbers;
    long long numComparisons;
    long long numSwaps;
    bool done;
};

// ---------------------------
// Bubble Sort Implementation
// ---------------------------
class BubbleSort : public SortAlgorithm {
public:
    BubbleSort(vector<int>& nums) 
        : SortAlgorithm(nums), passIndex(0), currentIndex(0) {}

    void step() override {
        if (passIndex >= static_cast<int>(numbers.size()) - 1) {
            done = true;
            return;
        }
        numComparisons++;
        if (numbers[currentIndex] > numbers[currentIndex + 1]) {
            swap(numbers[currentIndex], numbers[currentIndex + 1]);
            numSwaps++;
        }
        currentIndex++;
        if (currentIndex >= static_cast<int>(numbers.size()) - passIndex - 1) {
            currentIndex = 0;
            passIndex++;
        }
        if (isSorted(numbers))
            done = true;
    }

    bool isIndexActive(int idx) const override {
        return (idx == currentIndex || idx == currentIndex + 1);
    }

private:
    int passIndex, currentIndex;
};

// -----------------------------
// Insertion Sort Implementation
// -----------------------------
class InsertionSort : public SortAlgorithm {
public:
    InsertionSort(vector<int>& nums) 
        : SortAlgorithm(nums), i(1), j(0), key(0), holdingKey(false) {}

    void step() override {
        if (i >= static_cast<int>(numbers.size())) {
            done = true;
            return;
        }
        if (!holdingKey) {
            key = numbers[i];
            j = i - 1;
            holdingKey = true;
        }
        if (j >= 0)
            numComparisons++;
        if (j >= 0 && numbers[j] > key) {
            numbers[j + 1] = numbers[j];
            j--;
            numSwaps++;
        } else {
            numbers[j + 1] = key;
            i++;
            holdingKey = false;
            numSwaps++;
        }
    }

    bool isIndexActive(int idx) const override {
        return (idx == i || idx == j);
    }

private:
    int i, j, key;
    bool holdingKey;
};

// -------------------------
// Quick Sort Implementation
// -------------------------
struct Range {
    int left, right;
};

class QuickSort : public SortAlgorithm {
public:
    QuickSort(vector<int>& nums) 
        : SortAlgorithm(nums), inPartition(false), leftBound(0), rightBound(0), pivotIdx(0)
    {
        ranges.push_back({0, static_cast<int>(numbers.size()) - 1});
    }

    void step() override {
        if (!inPartition) {
            if (ranges.empty()) {
                done = true;
                return;
            }
            Range current = ranges.back();
            ranges.pop_back();
            leftBound = current.left;
            rightBound = current.right;
            pivotIdx = rightBound;
            inPartition = true;
            return;
        }
        int p = partition(leftBound, rightBound);
        if (p - 1 > leftBound)
            ranges.push_back({leftBound, p - 1});
        if (p + 1 < rightBound)
            ranges.push_back({p + 1, rightBound});
        inPartition = false;
        if (isSorted(numbers))
            done = true;
    }

    bool isIndexActive(int idx) const override {
        return (idx == pivotIdx);
    }

private:
    vector<Range> ranges;
    bool inPartition;
    int leftBound, rightBound;
    int pivotIdx;

    int partition(int left, int right) {
        int pivotVal = numbers[right];
        int i = left - 1;
        for (int j = left; j < right; j++) {
            numComparisons++;
            if (numbers[j] <= pivotVal) {
                i++;
                swap(numbers[i], numbers[j]);
                numSwaps++;
            }
        }
        swap(numbers[i + 1], numbers[right]);
        numSwaps++;
        return i + 1;
    }
};

// -------------------------
// Merge Sort Implementation
// -------------------------
class MergeSort : public SortAlgorithm {
public:
    MergeSort(vector<int>& nums) 
        : SortAlgorithm(nums), currentSize(1), merging(false), leftStart(0), mid(0), rightEnd(0)
    {
        temp.resize(numbers.size());
    }

    void step() override {
        if (currentSize >= static_cast<int>(numbers.size())) {
            done = true;
            return;
        }
        if (!merging) {
            if (leftStart >= static_cast<int>(numbers.size()) - 1) {
                currentSize *= 2;
                leftStart = 0;
                if (currentSize >= static_cast<int>(numbers.size()))
                    done = true;
                return;
            }
            merging = true;
            mid = min(leftStart + currentSize - 1, static_cast<int>(numbers.size()) - 1);
            rightEnd = min(leftStart + 2 * currentSize - 1, static_cast<int>(numbers.size()) - 1);
            return;
        } else {
            mergeSection(leftStart, mid, rightEnd);
            merging = false;
            leftStart += 2 * currentSize;
            if (isSorted(numbers))
                done = true;
        }
    }

    bool isIndexActive(int idx) const override {
        if (merging && idx >= leftStart && idx <= rightEnd)
            return true;
        return false;
    }

private:
    int currentSize;
    bool merging;
    int leftStart, mid, rightEnd;
    vector<int> temp;

    void mergeSection(int left, int mid, int right) {
        int i = left, j = mid + 1, k = left;
        while (i <= mid && j <= right) {
            numComparisons++;
            if (numbers[i] <= numbers[j])
                temp[k++] = numbers[i++];
            else {
                temp[k++] = numbers[j++];
                numSwaps++;
            }
        }
        while (i <= mid)
            temp[k++] = numbers[i++];
        while (j <= right)
            temp[k++] = numbers[j++];
        for (int idx = left; idx <= right; idx++)
            numbers[idx] = temp[idx];
    }
};

// -------------------------
// Heap Sort Implementation
// -------------------------
class HeapSort : public SortAlgorithm {
public:
    HeapSort(vector<int>& nums)
        : SortAlgorithm(nums), builtHeap(false), heapIndex(static_cast<int>(numbers.size()) / 2 - 1), heapSize(static_cast<int>(numbers.size()))
    {}

    void step() override {
        if (!builtHeap) {
            if (heapIndex < 0) {
                builtHeap = true;
                heapIndex = heapSize - 1;
                return;
            }
            heapify(heapSize, heapIndex);
            heapIndex--;
        } else {
            if (heapSize <= 1) {
                done = true;
                return;
            }
            swap(numbers[0], numbers[heapSize - 1]);
            numSwaps++;
            heapSize--;
            heapify(heapSize, 0);
        }
    }

    bool isIndexActive(int idx) const override {
        return (idx == 0);
    }

private:
    bool builtHeap;
    int heapIndex;
    int heapSize;

    void heapify(int n, int i) {
        int largest = i;
        int leftChild = 2 * i + 1;
        int rightChild = 2 * i + 2;
        if (leftChild < n) {
            numComparisons++;
            if (numbers[leftChild] > numbers[largest])
                largest = leftChild;
        }
        if (rightChild < n) {
            numComparisons++;
            if (numbers[rightChild] > numbers[largest])
                largest = rightChild;
        }
        if (largest != i) {
            swap(numbers[i], numbers[largest]);
            numSwaps++;
            heapify(n, largest);
        }
    }
};

// ---------------------------
// Selection Sort Implementation
// ---------------------------
class SelectionSort : public SortAlgorithm {
public:
    SelectionSort(vector<int>& nums)
        : SortAlgorithm(nums), i(0), j(0), minIndex(0)
    {}

    void step() override {
        if (i >= static_cast<int>(numbers.size()) - 1) {
            done = true;
            return;
        }
        if (j < static_cast<int>(numbers.size())) {
            numComparisons++;
            if (numbers[j] < numbers[minIndex])
                minIndex = j;
            j++;
            if (j >= static_cast<int>(numbers.size())) {
                if (minIndex != i) {
                    swap(numbers[i], numbers[minIndex]);
                    numSwaps++;
                }
                i++;
                j = i;
                minIndex = i;
            }
        }
    }

    bool isIndexActive(int idx) const override {
        return (idx == minIndex || idx == j);
    }

private:
    int i, j, minIndex;
};

// -------------------------
// Shell Sort Implementation
// -------------------------
class ShellSort : public SortAlgorithm {
public:
    ShellSort(vector<int>& nums)
        : SortAlgorithm(nums), gap(static_cast<int>(numbers.size()) / 2), i(gap), inLoop(false), currentJ(0), tempVal(0)
    {}

    void step() override {
        if (gap < 1) {
            done = true;
            return;
        }
        if (!inLoop) {
            if (i < static_cast<int>(numbers.size())) {
                tempVal = numbers[i];
                currentJ = i;
                inLoop = true;
            } else {
                gap = gap / 2;
                if (gap < 1) {
                    done = true;
                    return;
                }
                i = gap;
                return;
            }
        }
        if (inLoop) {
            if (currentJ >= gap) {
                numComparisons++;
                if (numbers[currentJ - gap] > tempVal) {
                    numbers[currentJ] = numbers[currentJ - gap];
                    numSwaps++;
                    currentJ -= gap;
                    return;
                }
            }
            numbers[currentJ] = tempVal;
            inLoop = false;
            i++;
        }
    }

    bool isIndexActive(int idx) const override {
        if (!inLoop) {
            return (idx == i);
        } else {
            return (idx == currentJ || (currentJ - gap >= 0 && idx == currentJ - gap));
        }
    }

private:
    int gap, i, currentJ, tempVal;
    bool inLoop;
};

// ---------------------------
// Cocktail Sort Implementation
// ---------------------------
class CocktailSort : public SortAlgorithm {
public:
    CocktailSort(vector<int>& nums)
        : SortAlgorithm(nums), start(0), end(static_cast<int>(numbers.size()) - 1), currentIdx(0), forward(true)
    {
        currentIdx = start;
    }

    void step() override {
        if (start >= end) {
            done = true;
            return;
        }
        if (forward) {
            if (currentIdx < end) {
                numComparisons++;
                if (numbers[currentIdx] > numbers[currentIdx + 1]) {
                    swap(numbers[currentIdx], numbers[currentIdx + 1]);
                    numSwaps++;
                }
                currentIdx++;
            } else {
                forward = false;
                end--;
                currentIdx = end;
            }
        } else {
            if (currentIdx > start) {
                numComparisons++;
                if (numbers[currentIdx - 1] > numbers[currentIdx]) {
                    swap(numbers[currentIdx - 1], numbers[currentIdx]);
                    numSwaps++;
                }
                currentIdx--;
            } else {
                forward = true;
                start++;
                currentIdx = start;
            }
        }
        if (isSorted(numbers))
            done = true;
    }

    bool isIndexActive(int idx) const override {
        if (forward)
            return (idx == currentIdx || idx == currentIdx + 1);
        else
            return (idx == currentIdx || idx == currentIdx - 1);
    }

private:
    int start, end, currentIdx;
    bool forward;
};

// ------------------------
// Comb Sort Implementation
// ------------------------
class CombSort : public SortAlgorithm {
public:
    CombSort(vector<int>& nums)
        : SortAlgorithm(nums), gap(static_cast<int>(numbers.size())), currentIdx(0), swapped(false)
    {}

    void step() override {
        if (currentIdx < static_cast<int>(numbers.size()) - gap) {
            numComparisons++;
            if (numbers[currentIdx] > numbers[currentIdx + gap]) {
                swap(numbers[currentIdx], numbers[currentIdx + gap]);
                numSwaps++;
                swapped = true;
            }
            currentIdx++;
        } else {
            if (gap > 1) {
                gap = static_cast<int>(gap / 1.3f);
                if (gap < 1)
                    gap = 1;
            }
            if (gap == 1 && !swapped)
                done = true;
            currentIdx = 0;
            swapped = false;
        }
    }

    bool isIndexActive(int idx) const override {
        return (idx == currentIdx || idx == currentIdx + gap);
    }

private:
    int gap, currentIdx;
    bool swapped;
};

// --------------------------
// Odd-Even Sort Implementation
// --------------------------
class OddEvenSort : public SortAlgorithm {
public:
    OddEvenSort(vector<int>& nums)
        : SortAlgorithm(nums), phase(0), offset(0), swapped(false)
    {}

    void step() override {
        if (isSorted(numbers)) {
            done = true;
            return;
        }
        int start = (phase == 0) ? 0 : 1;
        if (offset < static_cast<int>(numbers.size()) - 1) {
            int idx = start + offset;
            if (idx + 1 < static_cast<int>(numbers.size())) {
                numComparisons++;
                if (numbers[idx] > numbers[idx + 1]) {
                    swap(numbers[idx], numbers[idx + 1]);
                    numSwaps++;
                    swapped = true;
                }
            }
            offset++;
        } else {
            if (!swapped)
                done = true;
            offset = 0;
            phase = 1 - phase;
            swapped = false;
        }
    }

    bool isIndexActive(int idx) const override {
        int start = (phase == 0) ? 0 : 1;
        int current = start + offset;
        return (idx == current || idx == current + 1);
    }

private:
    int phase;
    int offset;
    bool swapped;
};

// ------------------------------------
// Class managing the visualization
// ------------------------------------
class SortingVisualizer {
public:
    SortingVisualizer()
        : window(nullptr),
          numBars(kDefaultBarCount),
          sortStepDelay(kDefaultSortDelay),
          isSorting(false),
          isPaused(false),
          sortingComplete(false),
          timerStopped(false),
          finalSortTime(0.0),
          currentAlgoIndex(0),
          finalAnimationIndex(0),
          lastFrameTimestamp(0.0f),
          waveMode(false)
    {
        numbers.resize(numBars);
        barHeights.resize(numBars, 0.0f);
    }

    ~SortingVisualizer() {
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glDeleteProgram(shaderProgram);
        glDeleteVertexArrays(1, &barVAO);
        glDeleteBuffers(1, &barVBO);
        glfwTerminate();
    }

    bool init() {
        if (!glfwInit()) {
            cerr << "Could not initialize GLFW\n";
            return false;
        }
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        window = glfwCreateWindow(kWindowWidth, kWindowHeight, "Sorting Visualizer", nullptr, nullptr);
        if (!window) {
            cerr << "Could not create GLFW window\n";
            glfwTerminate();
            return false;
        }
        glfwMakeContextCurrent(window);
        glfwSetFramebufferSizeCallback(window, framebufferResize);

        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
            cerr << "Failed to initialize GLAD\n";
            return false;
        }

        // Setup ImGui
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        setTheme();
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 330 core");

        setupShaders();
        setupBarGeometry();

        randomizeNumbers();
        resetAlgorithm();

        lastFrameTimestamp = static_cast<float>(glfwGetTime());
        return true;
    }

    void run() {
        bool spaceWasPressed = false;
        sortStartTime = chrono::steady_clock::now();
        lastStepTime = sortStartTime;
        finalAnimLastUpdate = sortStartTime;

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();

            bool spacePressed = (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS);
            if (spacePressed && !spaceWasPressed && isSorting)
                isPaused = !isPaused;
            spaceWasPressed = spacePressed;

            if (isSorting && !isPaused && !sortingComplete && !waveMode) {
                auto now = chrono::steady_clock::now();
                if (chrono::duration_cast<chrono::milliseconds>(now - lastStepTime) >= sortStepDelay) {
                    if (currentAlgorithm)
                        currentAlgorithm->step();
                    lastStepTime = now;
                    if (currentAlgorithm && currentAlgorithm->isDone() && !timerStopped) {
                        finalSortTime = chrono::duration_cast<chrono::milliseconds>(now - sortStartTime).count() / 1000.0;
                        timerStopped = true;
                        sortingComplete = true;
                    }
                }
            }

            updateBarHeights();

            if (sortingComplete && finalAnimationIndex < numBars) {
                auto now = chrono::steady_clock::now();
                if (chrono::duration_cast<chrono::milliseconds>(now - finalAnimLastUpdate).count() > 20) {
                    finalAnimationIndex++;
                    finalAnimLastUpdate = now;
                }
            }

            render();
            glfwSwapBuffers(window);
        }
    }

private:
    GLFWwindow* window;
    int numBars;
    vector<int> numbers;
    vector<float> barHeights;
    unique_ptr<SortAlgorithm> currentAlgorithm;
    chrono::steady_clock::time_point sortStartTime, lastStepTime, finalAnimLastUpdate;
    chrono::milliseconds sortStepDelay;
    bool isSorting, isPaused, sortingComplete, timerStopped;
    double finalSortTime;
    int currentAlgoIndex;
    int finalAnimationIndex;
    float lastFrameTimestamp;
    bool waveMode;
    vector<int> comboSequence;

    unsigned int shaderProgram, barVAO, barVBO;

    // -------------------------------
    // Data generation functions
    // -------------------------------
    void randomizeNumbers() {
        numbers.resize(numBars);
        barHeights.resize(numBars, 0.0f);
        for (int i = 0; i < numBars; i++) {
            numbers[i] = i + 1;
        }
        random_device rd;
        mt19937 gen(rd());
        shuffle(numbers.begin(), numbers.end(), gen);
    }

    void mountainNumbers() {
        numbers.resize(numBars);
        barHeights.resize(numBars, 0.0f);
        int mid = numBars / 2;
        for (int i = 0; i < numBars; i++) {
            if (i <= mid)
                numbers[i] = i + 1;
            else
                numbers[i] = numbers[mid - (i - mid)];
        }
        int currentMin = *min_element(numbers.begin(), numbers.end());
        int currentMax = *max_element(numbers.begin(), numbers.end());
        for (int i = 0; i < numBars; i++) {
            numbers[i] = 1 + (numbers[i] - currentMin) * (numBars - 1) / (currentMax - currentMin);
        }
    }

    void valleyNumbers() {
        numbers.resize(numBars);
        barHeights.resize(numBars, 0.0f);
        int mid = numBars / 2;
        for (int i = 0; i < numBars; i++) {
            if (i <= mid)
                numbers[i] = mid - i + 1;
            else
                numbers[i] = i - mid + 1;
        }
        int currentMin = *min_element(numbers.begin(), numbers.end());
        int currentMax = *max_element(numbers.begin(), numbers.end());
        for (int i = 0; i < numBars; i++) {
            numbers[i] = 1 + (numbers[i] - currentMin) * (numBars - 1) / (currentMax - currentMin);
        }
    }

    void reverseNumbers() {
        numbers.resize(numBars);
        barHeights.resize(numBars, 0.0f);
        for (int i = 0; i < numBars; i++) {
            numbers[i] = numBars - i;
        }
    }

    void resetAlgorithm() {
        sortingComplete = false;
        timerStopped = false;
        finalSortTime = 0.0;
        isPaused = false;
        finalAnimationIndex = 0;
        sortStartTime = chrono::steady_clock::now();
        lastStepTime = sortStartTime;
        finalAnimLastUpdate = sortStartTime;
        createAlgorithm(currentAlgoIndex);
    }

    void createAlgorithm(int index) {
        switch(index) {
            case 0: currentAlgorithm = make_unique<BubbleSort>(numbers); break;
            case 1: currentAlgorithm = make_unique<InsertionSort>(numbers); break;
            case 2: currentAlgorithm = make_unique<QuickSort>(numbers); break;
            case 3: currentAlgorithm = make_unique<MergeSort>(numbers); break;
            case 4: currentAlgorithm = make_unique<HeapSort>(numbers); break;
            case 5: currentAlgorithm = make_unique<SelectionSort>(numbers); break;
            case 6: currentAlgorithm = make_unique<ShellSort>(numbers); break;
            case 7: currentAlgorithm = make_unique<CocktailSort>(numbers); break;
            case 8: currentAlgorithm = make_unique<CombSort>(numbers); break;
            case 9: currentAlgorithm = make_unique<OddEvenSort>(numbers); break;
            default: currentAlgorithm = make_unique<BubbleSort>(numbers); break;
        }
    }

    void recordCombo(int code) {
        comboSequence.push_back(code);
        if(comboSequence.size() > 6)
            comboSequence.erase(comboSequence.begin());
        if(comboSequence.size() == 6 &&
           comboSequence[0] == 0 &&
           comboSequence[1] == 1 &&
           comboSequence[2] == 1 &&
           comboSequence[3] == 2 &&
           comboSequence[4] == 3 &&
           comboSequence[5] == 2) {
            waveMode = true;
            comboSequence.clear();
        }
    }

    // -------------------------------
    // OpenGL and ImGui setup functions
    // -------------------------------
    static void framebufferResize(GLFWwindow* window, int width, int height) {
        glViewport(0, 0, width, height);
    }

    unsigned int compileShader(unsigned int type, const char* source) {
        unsigned int shader = glCreateShader(type);
        glShaderSource(shader, 1, &source, nullptr);
        glCompileShader(shader);
        int success;
        char infoLog[512];
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if(!success) {
            glGetShaderInfoLog(shader, 512, nullptr, infoLog);
            cerr << "Shader compilation failed:\n" << infoLog << endl;
        }
        return shader;
    }

    void setupShaders() {
        const char* vertexSource = R"(
        #version 330 core
        layout(location=0) in vec2 aPos;
        uniform mat4 uMVP;
        void main(){
            gl_Position = uMVP * vec4(aPos, 0.0, 1.0);
        }
        )";
        const char* fragmentSource = R"(
        #version 330 core
        out vec4 FragColor;
        uniform vec4 uColor;
        void main(){
            FragColor = uColor;
        }
        )";
        unsigned int vertShader = compileShader(GL_VERTEX_SHADER, vertexSource);
        unsigned int fragShader = compileShader(GL_FRAGMENT_SHADER, fragmentSource);
        shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertShader);
        glAttachShader(shaderProgram, fragShader);
        glLinkProgram(shaderProgram);
        int success;
        char infoLog[512];
        glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
        if(!success) {
            glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog);
            cerr << "Shader program linking failed:\n" << infoLog << endl;
        }
        glDeleteShader(vertShader);
        glDeleteShader(fragShader);
    }

    void setupBarGeometry() {
        float vertices[] = {
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
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glBindVertexArray(0);
    }

    void setTheme() {
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

    // -------------------------------
    // Update bar heights for animation
    // -------------------------------
    void updateBarHeights() {
        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        float minHeight = 10.0f;
        float maxHeight = fbHeight * 0.9f;
        float currentTime = static_cast<float>(glfwGetTime());

        if(waveMode) {
            float base = (minHeight + maxHeight) / 2.0f;
            float amplitude = (maxHeight - minHeight) / 4.0f;
            float speed = 2.0f;
            float phaseOffset = 0.5f;
            for (int i = 0; i < numBars; i++) {
                barHeights[i] = base + amplitude * sin(currentTime * speed + i * phaseOffset);
            }
            return;
        }

        float dt = currentTime - lastFrameTimestamp;
        lastFrameTimestamp = currentTime;
        float interpFactor = 1 - exp(-10.0f * dt);

        for (int i = 0; i < numBars; i++) {
            float target = minHeight + (maxHeight - minHeight) * ((numbers[i] - 1) / float(numBars - 1));
            barHeights[i] = glm::mix(barHeights[i], target, interpFactor);
        }
    }

    // -------------------------------
    // Render the UI and bars
    // -------------------------------
    void render() {
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Sorting Visualizer");
        if (ImGui::BeginTabBar("Tabs")) {
            if (ImGui::BeginTabItem("Controls")) {
                if(waveMode) {
                    ImGui::Text("Wave Mode Active!");
                    if(ImGui::Button("Exit Wave Mode")) {
                        waveMode = false;
                        comboSequence.clear();
                        randomizeNumbers();
                        resetAlgorithm();
                        isSorting = false;
                    }
                } else if (!isSorting) {
                    if (ImGui::BeginCombo("Algorithm", algorithmNames[currentAlgoIndex])) {
                        for (int n = 0; n < IM_ARRAYSIZE(algorithmNames); n++) {
                            bool isSelected = (currentAlgoIndex == n);
                            if (ImGui::Selectable(algorithmNames[n], isSelected))
                                currentAlgoIndex = n;
                            if (isSelected)
                                ImGui::SetItemDefaultFocus();
                        }
                        ImGui::EndCombo();
                    }
                    if (ImGui::Button("Randomize")) {
                        randomizeNumbers();
                        resetAlgorithm();
                        recordCombo(0);
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Mountain")) {
                        mountainNumbers();
                        resetAlgorithm();
                        recordCombo(1);
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Valley")) {
                        valleyNumbers();
                        resetAlgorithm();
                        recordCombo(2);
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Reverse")) {
                        reverseNumbers();
                        resetAlgorithm();
                        recordCombo(3);
                    }
                    if (ImGui::Button("Start")) {
                        createAlgorithm(currentAlgoIndex);
                        resetAlgorithm();
                        isSorting = true;
                        isPaused = false;
                        sortStartTime = chrono::steady_clock::now();
                        lastStepTime = sortStartTime;
                    }
                } else {
                    if (!isPaused) {
                        if (ImGui::Button("Pause"))
                            isPaused = true;
                    } else {
                        if (ImGui::Button("Resume")) {
                            isPaused = false;
                            lastStepTime = chrono::steady_clock::now();
                        }
                        ImGui::SameLine();
                        if (ImGui::Button("Step"))
                            if (currentAlgorithm) currentAlgorithm->step();
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Reset")) {
                        resetAlgorithm();
                        isSorting = false;
                        isPaused = false;
                    }
                }
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Choose a pattern and start sorting");
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Settings")) {
                static float delayMs = float(sortStepDelay.count());
                if (ImGui::SliderFloat("Sort Speed (ms)", &delayMs, 1.0f, 100.0f, "%.0f"))
                    sortStepDelay = chrono::milliseconds(int(delayMs));
                static int bars = numBars;
                if (ImGui::SliderInt("Number of Bars", &bars, 10, 300)) {
                    if (bars != numBars) {
                        numBars = bars;
                        randomizeNumbers();
                        resetAlgorithm();
                    }
                }
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Statistics")) {
                ImGui::Text("Algorithm: %s", algorithmNames[currentAlgoIndex]);
                double elapsed = 0.0;
                if (sortingComplete && timerStopped)
                    elapsed = finalSortTime;
                else if (isSorting && !isPaused) {
                    auto now = chrono::steady_clock::now();
                    elapsed = chrono::duration_cast<chrono::milliseconds>(now - sortStartTime).count() / 1000.0;
                }
                ImGui::Text("Time: %.2f s", elapsed);
                if (currentAlgorithm) {
                    ImGui::Text("Comparisons: %lld", currentAlgorithm->getComparisons());
                    ImGui::Text("Swaps:       %lld", currentAlgorithm->getSwaps());
                }
                float progress = 0.0f;
                ImGui::ProgressBar(progress, ImVec2(-1.0f, 0.0f));
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("About")) {
                ImGui::BeginChild("AboutChild", ImVec2(400, 300), true, ImGuiWindowFlags_HorizontalScrollbar);
                ImGui::TextWrapped("%s", algorithmDescriptions[currentAlgoIndex]);
                ImGui::Spacing();
                ImGui::TextWrapped("Time Complexity:\nBubble Sort: O(n^2)\nInsertion Sort: O(n^2) (O(n) best-case)\nQuick Sort: O(n log n) average\nMerge Sort: O(n log n)\nHeap Sort: O(n log n)\nSelection Sort: O(n^2)\nShell Sort: ~O(n log^2 n)\nCocktail Sort: O(n^2)\nComb Sort: O(n^2) worst-case\nOdd-Even Sort: O(n^2)");
                ImGui::EndChild();
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
        ImGui::End();

        glClearColor(0.15f, 0.15f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        glm::mat4 projection = glm::ortho(0.0f, float(fbWidth), 0.0f, float(fbHeight), -1.0f, 1.0f);
        glUseProgram(shaderProgram);
        glBindVertexArray(barVAO);
        int locMVP = glGetUniformLocation(shaderProgram, "uMVP");
        int locColor = glGetUniformLocation(shaderProgram, "uColor");

        double mouseX, mouseY;
        glfwGetCursorPos(window, &mouseX, &mouseY);
        float oglMouseY = float(fbHeight) - float(mouseY);
        float barWidth = (fbWidth - (numBars + 1) * kBarGap) / float(numBars);

        for (int i = 0; i < numBars; i++) {
            float xPos = kBarGap + i * (barWidth + kBarGap);
            float height = barHeights[i];

            {
                glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(xPos + 2.0f, 2.0f, 0.0f));
                model = glm::scale(model, glm::vec3(barWidth, height, 1.0f));
                glm::mat4 mvp = projection * model;
                glUniformMatrix4fv(locMVP, 1, GL_FALSE, glm::value_ptr(mvp));
                glUniform4f(locColor, 0.1f, 0.1f, 0.1f, 0.5f);
                glDrawArrays(GL_TRIANGLES, 0, 6);
            }

            glm::vec4 color(0.9f, 0.9f, 0.9f, 1.0f);

            if (!sortingComplete) {
                if (currentAlgorithm && currentAlgorithm->isIndexActive(i)) {
                    color = glm::vec4(0.8f, 0.2f, 0.2f, 1.0f);
                }
            } else {
                color = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
                if (i <= finalAnimationIndex) {
                    float t = float(glfwGetTime()) * 5.0f;
                    float pulse = 0.5f + 0.5f * sin(t);
                    glm::vec4 base(0.4f, 0.9f, 0.4f, 1.0f);
                    float brightness = 0.8f + 0.2f * pulse;
                    color = glm::vec4(base.r * brightness, base.g * brightness, base.b * brightness, 1.0f);
                }
            }

            bool hover = (mouseX >= xPos && mouseX <= xPos + barWidth && oglMouseY <= height);
            if (hover) {
                color = glm::mix(color, glm::vec4(1.0f, 1.0f, 0.7f, 1.0f), 0.5f);
            }

            glUniform4fv(locColor, 1, &color[0]);

            glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(xPos, 0.0f, 0.0f));
            model = glm::scale(model, glm::vec3(barWidth, height, 1.0f));
            glm::mat4 mvp = projection * model;
            glUniformMatrix4fv(locMVP, 1, GL_FALSE, glm::value_ptr(mvp));
            glDrawArrays(GL_TRIANGLES, 0, 6);
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }
};

int main() {
    SortingVisualizer visualizer;
    if (!visualizer.init())
        return -1;
    visualizer.run();
    return 0;
}
