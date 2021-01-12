VERSION 5.00
Begin VB.Form Form1 
   Caption         =   """fbench"" Floating-point benchmark"
   ClientHeight    =   4215
   ClientLeft      =   60
   ClientTop       =   510
   ClientWidth     =   8010
   Icon            =   "fbench.frx":0000
   LinkTopic       =   "Form1"
   ScaleHeight     =   4215
   ScaleWidth      =   8010
   StartUpPosition =   3  'Windows Default
   Begin VB.CommandButton Command1 
      Caption         =   "Run"
      BeginProperty Font 
         Name            =   "Comic Sans MS"
         Size            =   9.75
         Charset         =   0
         Weight          =   700
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   375
      Left            =   60
      TabIndex        =   3
      Top             =   45
      Width           =   840
   End
   Begin VB.TextBox Text1 
      Alignment       =   2  'Center
      BeginProperty Font 
         Name            =   "Comic Sans MS"
         Size            =   9.75
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   375
      Left            =   6210
      TabIndex        =   1
      Text            =   "1000"
      Top             =   75
      Width           =   1005
   End
   Begin VB.ListBox List1 
      BackColor       =   &H00FFFFF8&
      BeginProperty Font 
         Name            =   "Lucida Console"
         Size            =   12
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   1740
      Left            =   4110
      TabIndex        =   0
      Top             =   495
      Width           =   2850
   End
   Begin VB.Label Label1 
      Alignment       =   1  'Right Justify
      Caption         =   "Iterations"
      BeginProperty Font 
         Name            =   "Comic Sans MS"
         Size            =   9
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   300
      Left            =   4755
      TabIndex        =   2
      Top             =   135
      Width           =   1395
   End
End
Attribute VB_Name = "Form1"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit

'===================================================================
'  John Walker's FBENCH benchmark test
'
'  VB6 adaptation by Jim White, mathimagics@yahoo.co.uk
'
'  Operation:
'       adjust the iteration count in the TextBox as desired
'       and press the RUN button
'
'  Compiler switches:
'
'       fbench.exe was compiled with the VB6 default options.
'       fbench_opt.exe was made with the default "runtime
'       checks" disabled, viz.:
'
'            - Array bounds
'            - Floating Point exceptions
'            - Integer overflow
'
'
'  Typical results:
'       Times for 1000 iterations on a Pentium 4 (3.0Ghz)
'         based on an actual run of 1,000,000 iterations
'
'       EXE with default options:  (Compile for "Fast execution")
'
'            0.01852 secs
'
'       EXE with Fast Execution + Remove array bounds checks + Remove
'            Integr Overflow Checks + Remove Floating Point Checks
'
'            0.01307 secs
'
'       EXE with INTRIG = False  (using inline trig functions)
'
'            0.02975 secs
'
'       C version (gcc with -O3):
'
'            0.01435 secs
'
'       IDE mode:
'            0.0859  secs
'
'===================================================================
'
'  NOTE particularly which version was the fastest overall?  :)
'
'  Next time you encounter the popular mythology regarding how
'  terribly slow VB6 is compared to a "real language" (!) like
'  C, you can show them the errors in their thinking.
'
'  The VBP file I have supplied has the default compiler
'  options specified. Go to menu option
'        Project Properties -> Compile -> Advanced Optimisations
'  in order to access the "Remove Checks" compiler switches
'
'  Jim White
'  mathimagics@yahoo.co.uk
'  28 July 2005
'
'===================================================================

Private Sub Form_Resize()
   If WindowState = 1 Then Exit Sub
   List1.Move 30, List1.Top, Me.ScaleWidth - 60, Me.ScaleHeight - List1.Top - 60
   Text1.Left = Me.ScaleWidth - Text1.Width - 300
   Label1.Left = Text1.Left - Label1.Width - 60
   End Sub

Private Sub Command1_Click()
   fbench_run
   End Sub

