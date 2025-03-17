import numpy as np
from scipy.special import factorial, lpmn
import vispy.scene
from PyQt5 import QtWidgets, QtCore

def Laguerre(alpha, degree, x): 
    if degree == 0:
        return 1 
    if degree == 1:
        return -x + alpha + 1   
    else:
        L = 0
        L0, L1 = 1, -x + alpha + 1
        
        for i in range(2, degree + 1):   
            L = (2 + (alpha - 1 - x) / i) * L1 - (1 + (alpha - 1) / i) * L0
            L0, L1 = L1, L
        return L
    
def Legendre(m,l,x):
    return lpmn(m,l,x)[0][m][l]   
Legendre_ = np.vectorize(Legendre, excluded=['m', 'l'])


def Radial(r,n,l):
    rho = 2 * r / n
    return (2 / n**2) * np.sqrt(factorial(n - l - 1) / factorial(n + l)) * (rho**l) * Laguerre(2 * l + 1, n - l - 1, rho) * np.exp(-rho / 2)

def Ylm(t,l,m):
    if m >= 0:
        N = np.sqrt((2 * l + 1) * factorial((l - m)) / (4 * np.pi * factorial((l + m))))
        return N * Legendre_(m, l, np.cos(t)) * (-1)**m
    else:
        N = np.sqrt((2 * l + 1) * factorial((l + m)) / (4 * np.pi * factorial((l - m))))
        return N * Legendre_(-m, l, np.cos(t))  


length = 50
pixels = 100
x = np.linspace(-length, length, pixels)
y = np.linspace(-length, length, pixels)
z = np.linspace(-length, length, pixels)

X, Y, Z = np.meshgrid(x, y, z)

R = np.sqrt(X**2 + Y**2 + Z**2)
theta = np.arccos(Z / R)
phi = np.arctan2(Y, X)

n = 5
l = 3
m = 1

F = (Radial(R, n, l) * Ylm(theta, l, m) * np.abs(np.exp(1j * m * phi)))**2  

max=100000
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Interactive Isosurface Level : (n, l, m) = ('+str(n)+', '+str(l)+', '+str(m)+')')
        self.canvas = vispy.scene.SceneCanvas(keys='interactive', show=True)
        self.view = self.canvas.central_widget.add_view()

        # Create an Isosurface object
        self.surface = vispy.scene.visuals.Isosurface(F, level=F.max()/10., color=(1, 1, 0, 1), shading='smooth', parent=self.view.scene)
        self.surface.transform = vispy.scene.transforms.STTransform(translate=(-pixels/2, -pixels/2, -pixels/2))

        # Configure the camera
        self.cam = vispy.scene.TurntableCamera(elevation=30, azimuth=30)
        self.cam.set_range((-10, 10), (-10, 10), (-10, 10))
        self.view.camera = self.cam

        # Add the canvas to the main window
        self.setCentralWidget(self.canvas.native)

        # Add a slider to adjust the level parameter
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider.setMinimum(1)
        self.slider.setMaximum(max)
        self.slider.valueChanged.connect(self.update_level)

        # Add the slider to the main window
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.canvas.native)
        self.layout.addWidget(self.slider)
        
        self.central_widget = QtWidgets.QWidget()
        self.central_widget.setLayout(self.layout)
        self.setCentralWidget(self.central_widget)

    def update_level(self):
        level_value = F.max() * (self.slider.value()/max)
        self.surface.level = level_value
        self.surface.update()  
        self.canvas.update()
        print(level_value)


def main():
    appQt = QtWidgets.QApplication([])
    win = MainWindow()
    win.show()
    appQt.exec_()

if __name__ == '__main__':
    main()
