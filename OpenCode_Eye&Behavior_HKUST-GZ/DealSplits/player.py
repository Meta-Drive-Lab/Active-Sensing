from backupFunctions import *
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QLabel, QFileDialog, QSlider, QHBoxLayout, QVBoxLayout, QWidget,QLineEdit
from PyQt5.QtCore import QTimer, Qt


"""
fixitem 1:road 2:mirror
"""

"""
the item value cases
road
     3: The car in front is very close and cannot change lanes
     2: There is no car ahead and you can change lanes
     1: The car in front is far away and you can change lanes

mirror
     3: There are no cars on the side and you can change lanes
     2: There is a car on the side, it is far away, and you can change lanes
     1: There is a car on the side. It is very close and cannot change the lane.
"""


'''
This is a video player code that allows you to display the gaze position and predicted rear view mirror position in the image
You can set the value of itemvalue according to the specific conditions of the road according to the several situations mentioned above.
You can give the value of choice by whether to change lanes
You can give a start frame and a stop frame and re-intercept the interval
For other small functions, please refer to the function implementation.
'''
classes = ['R', 'fan', 'armor', 'mid', 'left', 'right']

class VideoPlayer(QMainWindow):
    
    # Constructor for the VideoPlayer class
    # Initializes the user interface elements, variables, and sets up the video player.
    def __init__(self):
        super().__init__()

        # set up a window
        self.setWindowTitle('Video Player')

        # Create a Timer object
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.play)

        # Create OpenCV video capture object
        self.cap = None
        self.etgdatas = []
        self.video_files = sorted([f for f in os.listdir('./videos') if f.endswith(('.avi', '.mp4', '.flv', '.mkv'))])
        self.current_video_index = 0
        etgname = os.path.splitext(self.video_files[self.current_video_index])[0] + ".txt"
        self.etgpath = os.path.join('./etgdatas', etgname)
        self.outpath = os.path.join('./out', etgname)
        if os.path.exists(self.etgpath):
            print(self.etgpath)

        # Create buttons
        self.btn_open = QPushButton('Open Video', self)
        self.btn_open.clicked.connect(self.open_video)
        
        self.btn_play = QPushButton('Play', self)
        self.btn_play.clicked.connect(self.start_play)

        self.btn_pause = QPushButton('Pause', self)
        self.btn_pause.clicked.connect(self.pause_play)

        self.btn_pre = QPushButton('Previous Video', self)
        self.btn_pre.clicked.connect(self.pre_video)

        self.btn_next = QPushButton('Next Video', self)
        self.btn_next.clicked.connect(self.next_video)

        self.btn_delete = QPushButton('Delete Video', self)
        self.btn_delete.clicked.connect(self.delete_video)

        self.btn_save = QPushButton('Save', self)
        self.btn_save.clicked.connect(self.save_endframe)

        self.btn_iflanechange_yes = QPushButton('Yes', self)
        self.btn_iflanechange_yes.clicked.connect(self.write_iflanechange_yes)

        self.btn_iflanechange_no = QPushButton('No', self)
        self.btn_iflanechange_no.clicked.connect(self.write_iflanechange_no)

        self.road_input = QLineEdit(self)
        self.mirror_input = QLineEdit(self)


        # Create a label to display the frame number
        self.label = QLabel(self)
        # Create a label to display the current video file name
        self.video_name_label = QLabel(self)
        
        if self.cap != None:
            self.update_label()
            self.update_video_name_label()

        # Create a slider
        self.slider = QSlider(Qt.Horizontal, self)
        self.slider.sliderMoved.connect(self.set_position)

        # Set layout
        hlayout = QHBoxLayout()
        hlayout.addWidget(self.btn_open)
        hlayout.addWidget(self.btn_play)
        hlayout.addWidget(self.btn_pause)
        hlayout.addWidget(self.btn_pre)
        hlayout.addWidget(self.btn_next)
        hlayout.addWidget(self.label)
        
        hlayout.addWidget(self.video_name_label)

        hlayout2 = QHBoxLayout()
        hlayout2.addWidget(self.btn_delete)
        hlayout2.addStretch()  

        hlayout3 = QHBoxLayout()
        hlayout3.addWidget(self.btn_save)
        hlayout3.addWidget(QLabel("Lane Change: "))
        hlayout3.addWidget(self.btn_iflanechange_yes)
        hlayout3.addWidget(self.btn_iflanechange_no)
        hlayout3.addStretch()

        hlayout4 = QHBoxLayout()
        hlayout4.addWidget(QLabel("Road: "))
        hlayout4.addWidget(self.road_input)
        hlayout4.addWidget(QLabel("Mirror: "))
        hlayout4.addWidget(self.mirror_input)


        vlayout = QVBoxLayout()
        vlayout.addLayout(hlayout)
        vlayout.addLayout(hlayout2)
        vlayout.addLayout(hlayout3)
        vlayout.addLayout(hlayout4)
        vlayout.addWidget(self.slider)

        central_widget = QWidget()
        central_widget.setLayout(vlayout)
        self.setCentralWidget(central_widget)

        self.show()


    # Opens and initializes a video file.
    # Displays the gaze position and rearview mirror position.
    def open_video(self):
        videopath = os.path.join('./videos', self.video_files[self.current_video_index])
        self.cap = cv2.VideoCapture(videopath)
        etgname = os.path.splitext(self.video_files[self.current_video_index])[0] + ".txt"
        self.etgpath = os.path.join('./etgdatas', etgname)
        self.outpath = os.path.join('./out', etgname)
        if os.path.exists(self.etgpath):
            print(self.etgpath)
        with open(self.outpath, "r") as f:
            lines = f.readlines()
        choicedata = lines[1].strip()
        if choicedata == 'no':
            self.next_video()
        self.etgdatas = deallabels(self.etgpath)
        self.slider.setMaximum(int(self.cap.get(cv2.CAP_PROP_FRAME_COUNT) - 1))
        self.slider.setValue(0)
        self.update_video_name_label()


    # Plays the video frame by frame.
    # Updates the display to show the video frames and gaze positions.
    def play(self):
        if self.cap:
            ret, frame = self.cap.read()
            if ret:
                self.update_label()
                self.slider.setValue(self.frame_number)
                if self.frame_number in self.etgdatas:
                    for x, y in self.etgdatas[self.frame_number]:
                        cv2.circle(frame, (x, y), 5, (255, 0, 0), -1)
                cv2.imshow("Video", frame)
            else:
                self.next_video()  # The current video has finished playing and the next video is loaded.

    # Starts playing the video from the current frame.
    def start_play(self):
        if self.cap:
            self.timer.start(30)

    # Loads and plays the next video in the list.
    def next_video(self):
        self.next_videoplay()
        self.play()

    # Loads and plays the previous video in the list.
    def pre_video(self):
        self.pre_videoplay()
        self.play()

    # Pauses the video playback.
    def pause_play(self):
        self.timer.stop()

    # Sets the video playback position based on the slider's position.
    def set_position(self, position):
        if self.cap:
            self.cap.set(cv2.CAP_PROP_POS_FRAMES, position)
            
            self.update_label()
            ret, frame = self.cap.read()
            if ret:
                if self.frame_number in self.etgdatas:
                    for x, y in self.etgdatas[self.frame_number]:
                        cv2.circle(frame, (x, y), 5, (255, 0, 0), -1)
                cv2.imshow("Video", frame)

    # Updates the label to display the current frame number.
    def update_label(self):
        self.frame_number = int(self.cap.get(cv2.CAP_PROP_POS_FRAMES))
        self.label.setText(f"Frame: {self.frame_number}")

    # Updates the label to display the name of the current video.
    def update_video_name_label(self):
        current_video_name = self.video_files[self.current_video_index] if self.video_files else "No Video"
        self.video_name_label.setText(f"Video: {current_video_name}")

    # Stops the current video, goes to the previous video, and starts playing it.
    def pre_videoplay(self):
        self.cap.release()
        cv2.destroyAllWindows()
        self.current_video_index -= 1
        if self.current_video_index >= 0:
            self.open_video()
        else:
            self.current_video_index = 0
            self.open_video()

    # Saves the user's input for road and mirror values to the output file.
    def save_endframe(self):
        pass
        road_value = self.road_input.text()
        mirror_value = self.mirror_input.text()
        with open(self.outpath, 'a') as f:
            f.write(f"road: {road_value}\n")
            f.write(f"mirror: {mirror_value}\n")
        
    # Writes "yes" to indicate a lane change in the output file.
    def write_iflanechange_yes(self):
        # pass
        with open(self.outpath, 'a') as f:
            f.write("yes\n")

    # Writes "no" to indicate no lane change in the output file.        
    def write_iflanechange_no(self):
        # pass
        with open(self.outpath, 'a') as f:
            f.write("no\n")

    # Stops the current video, goes to the next video, and starts playing it.
    def next_videoplay(self):
        self.cap.release()
        cv2.destroyAllWindows()
        self.current_video_index += 1
        if self.current_video_index < len(self.video_files):
            self.open_video()
        else:
            self.current_video_index = 0
            self.open_video()


    # Deletes the current video and its associated files.
    # Handles updating the UI and moving to the next video.
    def delete_video(self):
        if self.cap:
            self.timer.stop()  
            self.cap.release()
            cv2.destroyAllWindows()
            video_file_path = os.path.join('./videos', self.video_files[self.current_video_index])
            etgname = os.path.splitext(self.video_files[self.current_video_index])[0] + ".txt"
            self.etgpath = os.path.join('./etgdatas', etgname)
            self.outpath = os.path.join('./out', etgname)
            if os.path.exists(self.etgpath):
                os.remove(self.etgpath)
                os.remove(self.outpath)
            os.remove(video_file_path)
            self.video_files.pop(self.current_video_index) 

            if not self.video_files:  
                self.cap = None
                self.slider.setValue(0)
                self.label.setText("Frame: 0")
                self.update_video_name_label()
            else:  
                if self.current_video_index == len(self.video_files):
                    self.current_video_index -= 1
                self.open_video()

# Entry point of the application
if __name__ == '__main__':
    app = QApplication(sys.argv)
    player = VideoPlayer()
    sys.exit(app.exec_())
