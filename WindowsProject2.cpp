#include <windows.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <unordered_map>

#define M_PI 3.14159

struct vec3 {
    float x, y, z;
    float dx, dy, dz;
    float g_dy;
    float density;
    float pressure;
};

// Global variables
HBITMAP hBitmap = NULL;
const int bitmapWidth = 600; // Width of the buffer
const int bitmapHeight = 600; // Height of the buffer
const int bitmapDepth = 600;

const int particleBoxSizeX = 20;
const int particleBoxSizeY = 20;
const int particleBoxSizeZ = 20;
const int numberOfBoxes = (bitmapHeight / particleBoxSizeY) * (bitmapWidth / particleBoxSizeX) * (bitmapWidth / particleBoxSizeZ);
const int groundY = bitmapHeight;

const float H = 50;

// Camera position and orientation
float cameraX = bitmapWidth / 2.0f;
float cameraY = -bitmapHeight * 2.0f; // Higher than the scene
float cameraZ = bitmapHeight / 0.5f;

// Rotation angle (45 degrees in radians)
float yAngle = 0.0f * M_PI / 180.0f;
float angle = 45.0f * M_PI / 180.0f;

std::vector<vec3> particles[numberOfBoxes];

std::vector<float> kernelTable;
std::vector<float> dWdrTable;

void generate_particles() {
    for (int i = 0; i < 1000; i++) {
        int r1 = bitmapWidth / 3 + rand() % bitmapWidth/30;
        int r2 = bitmapHeight / 3 + rand() % bitmapHeight/30;
        int r3 = bitmapDepth / 3 + rand() % bitmapDepth/30;
        vec3 vec = { static_cast<float>(r1), static_cast<float>(r2), static_cast<float>(r3) };
        int R1 = r1 / particleBoxSizeX;
        int R2 = r2 / particleBoxSizeY;
        int R3 = r3 / particleBoxSizeZ;
        int index = R1 + R2 * (bitmapWidth / particleBoxSizeX) + R3 * (bitmapWidth / particleBoxSizeZ);
        particles[index].push_back(vec);
    }
}

float r(vec3 p1, vec3 p2) {
    return pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2);
}

float W(float r, float h) {
    if (r < h / 2) {
        return ((8 / (M_PI * pow(h, 3))) * (1 - 6 * pow(r / h, 2) + 6 * pow(r / h, 3)));
    }
    else if (r < h) {
        return (8 / (M_PI * pow(h, 3))) * 2 * pow(1 - r / h, 3);
    }
    else {
        return 0;
    }
}

float dW_dr(float r, float h) {
    if (r < h / 2) {
        return (8 / (M_PI * pow(h, 4))) * (-12 * r / h + 18 * pow(r, 2) / pow(h, 2));
    }
    else if (r < h) {
        return (8 / (M_PI * pow(h, 4))) * (-6 * pow(1 - r / h, 2) / h);
    }
    else {
        return 0;
    }
}

void initializeKernelTables(float h, int resolution) {
    kernelTable.resize(resolution);
    dWdrTable.resize(resolution);

    float stepSize = h / resolution;
    for (int i = 0; i < resolution; ++i) {
        float r = i * stepSize;
        kernelTable[i] = W(r, h);
        dWdrTable[i] = dW_dr(r, h);
    }
}

float W_lookup(float r, float h) {
    int index = static_cast<int>(r / h * kernelTable.size());
    index = max(0, min(index, static_cast<int>(kernelTable.size()) - 1));
    return kernelTable[index];
}

float dWdr_lookup(float r, float h) {
    int index = static_cast<int>(r / h * dWdrTable.size());
    index = max(0, min(index, static_cast<int>(dWdrTable.size()) - 1));
    return dWdrTable[index];
}

vec3 gradW(vec3 p1, vec3 p2, float h) {
    vec3 r = { p1.x - p2.x, p1.y - p2.y, p1.z - p2.z };
    float r_magnitude = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
    if (r_magnitude == 0) return { 0, 0, 0 };
    float dWdr = dWdr_lookup(r_magnitude, h);
    return { dWdr * (r.x / r_magnitude), dWdr * (r.y / r_magnitude), dWdr * (r.z / r_magnitude) };
}

void update_density() {
    for (int j = 0; j < numberOfBoxes; j++) {
        for (size_t i = 0; i < particles[j].size(); ++i) {
            vec3& p = particles[j][i];
            p.density = 0;
            for (int offsetZ = -1; offsetZ <= 1; ++offsetZ) {
                for (int offsetY = -1; offsetY <= 1; ++offsetY) {
                    for (int offsetX = -1; offsetX <= 1; ++offsetX) {
                        int neighborX = (j % (bitmapWidth / particleBoxSizeX)) + offsetX;
                        int neighborY = ((j / (bitmapWidth / particleBoxSizeX)) % (bitmapHeight / particleBoxSizeY)) + offsetY;
                        int neighborZ = (j / (bitmapWidth / particleBoxSizeX) / (bitmapHeight / particleBoxSizeY)) + offsetZ;

                        if (neighborX >= 0 && neighborX < (bitmapWidth / particleBoxSizeX) &&
                            neighborY >= 0 && neighborY < (bitmapHeight / particleBoxSizeY) &&
                            neighborZ >= 0 && neighborZ < (bitmapWidth / particleBoxSizeZ)) {
                            int neighborIndex = neighborX + neighborY * (bitmapWidth / particleBoxSizeX) + neighborZ * (bitmapWidth / particleBoxSizeZ);

                            for (size_t k = 0; k < particles[neighborIndex].size(); ++k) {
                                const vec3& p2 = particles[neighborIndex][k];
                                p.density += W_lookup(r(p, p2), H);
                            }
                        }
                    }
                }
            }
        }
    }
}

void update_pressure() {
    for (int j = 0; j < numberOfBoxes; j++) {
        for (size_t i = 0; i < particles[j].size(); ++i) {
            vec3& p_i = particles[j][i];
            vec3 force = { 0, 0, 0 };
            for (int offsetZ = -1; offsetZ <= 1; ++offsetZ) {
                for (int offsetY = -1; offsetY <= 1; ++offsetY) {
                    for (int offsetX = -1; offsetX <= 1; ++offsetX) {
                        int neighborX = (j % (bitmapWidth / particleBoxSizeX)) + offsetX;
                        int neighborY = ((j / (bitmapWidth / particleBoxSizeX)) % (bitmapHeight / particleBoxSizeY)) + offsetY;
                        int neighborZ = (j / (bitmapWidth / particleBoxSizeX) / (bitmapHeight / particleBoxSizeY)) + offsetZ;

                        if (neighborX >= 0 && neighborX < (bitmapWidth / particleBoxSizeX) &&
                            neighborY >= 0 && neighborY < (bitmapHeight / particleBoxSizeY) &&
                            neighborZ >= 0 && neighborZ < (bitmapWidth / particleBoxSizeZ)) {
                            int neighborIndex = neighborX + neighborY * (bitmapWidth / particleBoxSizeX) + neighborZ * (bitmapWidth / particleBoxSizeZ);

                            for (size_t k = 0; k < particles[neighborIndex].size(); ++k) {
                                const vec3& p_j = particles[neighborIndex][k];
                                float pressure_i = p_i.density - 3;
                                float pressure_j = p_j.density - 3;
                                float density_j = p_j.density;

                                vec3 gradW_ij = gradW(p_i, p_j, H);
                                force.x += (pressure_i + pressure_j) / (2.0f * density_j) * gradW_ij.x;
                                force.y += (pressure_i + pressure_j) / (2.0f * density_j) * gradW_ij.y;
                                force.z += (pressure_i + pressure_j) / (2.0f * density_j) * gradW_ij.z;
                            }
                        }
                    }
                }
            }
            p_i.dx = 5 * force.x;
            p_i.dy = 5 * (force.y);
            p_i.dz = 5 * force.z;
        }
    }
}

void update_particles() {
    update_density();
    update_pressure();

    std::vector<std::vector<vec3>> updatedParticles(numberOfBoxes);

    for (int j = 0; j < numberOfBoxes; j++) {
        for (size_t i = 0; i < particles[j].size(); ++i) {
            vec3& p1 = particles[j][i];
            p1.x += p1.dx * 1;
            p1.y += p1.dy * 1;
            p1.z += p1.dz * 1;

            if (p1.x < 0) p1.x = 0;
            if (p1.x >= bitmapWidth) p1.x = bitmapWidth - 1;
            if (p1.y < 0) p1.y = 0;
            if (p1.y >= groundY) {
                p1.y = groundY - 1;
                p1.g_dy = 0;
            }
            if (p1.z < 0) p1.z = 0;
            if (p1.z >= bitmapDepth) p1.z = bitmapDepth - 1;

            int R1 = p1.x / particleBoxSizeX;
            int R2 = p1.y / particleBoxSizeY;
            int R3 = p1.z / particleBoxSizeZ;
            int index = R1 + R2 * (bitmapWidth / particleBoxSizeX) + R3 * (bitmapWidth / particleBoxSizeZ);
            updatedParticles[index].push_back(p1);
        }
    }
    for (int j = 0; j < numberOfBoxes; j++) {
        particles[j] = updatedParticles[j];
    }
}

void draw_buffer(HDC hdcBuffer) {
    // Create a BITMAPINFO structure to describe the bitmap format
    BITMAPINFO bmi = {};
    bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bmi.bmiHeader.biWidth = bitmapWidth;
    bmi.bmiHeader.biHeight = -bitmapHeight; // Negative height to indicate a top-down bitmap
    bmi.bmiHeader.biPlanes = 1;
    bmi.bmiHeader.biBitCount = 24; // 24-bit color depth (RGB)
    bmi.bmiHeader.biCompression = BI_RGB;

    // Create a buffer to hold pixel data (24 bits per pixel)
    BYTE* buffer = new BYTE[bitmapWidth * bitmapHeight * 3]; // 3 bytes per pixel (RGB)

    // Clear the buffer to white
    memset(buffer, 255, bitmapWidth * bitmapHeight * 3); // White color (RGB: 255, 255, 255)

    // Perspective projection parameters
    float fov = 90.0f; // Field of view in degrees
    float aspectRatio = static_cast<float>(bitmapWidth) / static_cast<float>(bitmapHeight);
    float nearPlane = 300.0f; // Distance to the nar plane
    float scale = tan(fov * 0.5 * M_PI / 180.0) * nearPlane;


    float cosAngle = cos(angle);
    float sinAngle = sin(angle);



    // Fixed depth range for z axis
    float minZ = 0.0f;
    float maxZ = bitmapHeight;

    // Draw particles onto the buffer
    for (int i = 0; i < numberOfBoxes; ++i) {
        for (size_t j = 0; j < particles[i].size(); ++j) {
            vec3& p = particles[i][j];

            // Translate particle position to camera space
            float relX = p.x - cameraX;
            float relY = p.y - cameraY;
            float relZ = p.z - cameraZ;

            // Apply rotation around X-axis to simulate a 45-degree view
            float yRotated = relY * cosAngle - relZ * sinAngle;
            float zRotated = relY * sinAngle + relZ * cosAngle;

            float xRotated = relX * cos(yAngle) + zRotated * sin(yAngle);
            float zRotatedY = -relX * sin(yAngle) + zRotated * cos(yAngle);

            // Perspective projection transformation
            if (yRotated != 0) { // Prevent division by zero
                float xProjected = (relX / yRotated) * scale * aspectRatio + bitmapWidth / 2.0f;
                float zProjected = (zRotated / yRotated) * scale + bitmapHeight / 2.0f;

                // Normalize zProjected to fit within the depth range
                float zNormalized = (p.z - minZ) / (maxZ - minZ);

                // Check if particle is within buffer bounds
                if (xProjected >= 0 && xProjected < bitmapWidth && zProjected >= 0 && zProjected < bitmapHeight) {
                    int offset = (static_cast<int>(zProjected) * bitmapWidth + static_cast<int>(xProjected)) * 3; // 3 bytes per pixel

                    // Set pixel color to black at particle coordinates
                    buffer[offset] = 0;       // Blue component
                    buffer[offset + 1] = 0;   // Green component
                    buffer[offset + 2] = static_cast<BYTE>(255 * zNormalized); // Red component (depth-based coloring)
                }
            }
        }
    }

    // Set the pixel data in the buffer to the device context, covering the entire client area
    SetDIBitsToDevice(hdcBuffer, 0, 0, bitmapWidth, bitmapHeight, 0, 0, 0, bitmapHeight, buffer, &bmi, DIB_RGB_COLORS);

    // Clean up
    delete[] buffer;
}

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
    switch (uMsg) {
    case WM_CREATE:
    {
        HDC hdc = GetDC(hwnd);
        hBitmap = CreateCompatibleBitmap(hdc, bitmapWidth, bitmapHeight);
        ReleaseDC(hwnd, hdc);
    }
    break;
    case WM_PAINT:
    {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);

        HDC hdcBuffer = CreateCompatibleDC(hdc);
        SelectObject(hdcBuffer, hBitmap);
        draw_buffer(hdcBuffer);
        BitBlt(hdc, 0, 0, bitmapWidth, bitmapHeight, hdcBuffer, 0, 0, SRCCOPY);
        DeleteDC(hdcBuffer);

        EndPaint(hwnd, &ps);
    }
    break;
    case WM_DESTROY:
        if (hBitmap != NULL)
            DeleteObject(hBitmap);
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hwnd, uMsg, wParam, lParam);
    }
    return 0;
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PSTR lpCmdLine, int nCmdShow) {
    AllocConsole();
    FILE* pConsoleIn, * pConsoleOut;
    freopen_s(&pConsoleIn, "CONIN$", "r", stdin);
    freopen_s(&pConsoleOut, "CONOUT$", "w", stdout);
    freopen_s(&pConsoleOut, "CONOUT$", "w", stderr);
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    std::cout << "hi1" << std::endl;
    generate_particles();
    std::cout << "hi1" << std::endl;

    initializeKernelTables(20, 1000);
    std::cout << "hi1" << std::endl;

    const wchar_t CLASS_NAME[] = L"SampleWindowClass";

    WNDCLASS wc = {};
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;

    RegisterClass(&wc);

    HWND hwnd = CreateWindowEx(
        0,
        CLASS_NAME,
        L"My Window",
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, bitmapWidth, bitmapHeight,
        NULL, NULL, hInstance, NULL
    );

    if (hwnd == NULL)
        return 0;

    ShowWindow(hwnd, nCmdShow);

    MSG msg = {};
    //Sleep(8000);
    while (true) {
        while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
            if (msg.message == WM_QUIT)
                return 0;
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
        std::cout << "hi" << std::endl;
        update_particles();
        InvalidateRect(hwnd, NULL, TRUE);
        UpdateWindow(hwnd);
    }

    return 0;
}
