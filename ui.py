import PySimpleGUI as sg
import os   
import yaml        as ya
from   pathlib     import Path

###############################################################################
#Non-standard Imports
###############################################################################
try:
    import SENAX as snx
except Exception as e:
    if Path(__file__).parent == Path.cwd():
        import addpath  
        import SENAX as snx
    else:
        raise e

###############################################################################
#HiPADGUI class
###############################################################################
class HiPADGUI(sg.Window):
    def __init__(self):

        #Menu
        menu_def = [['&Help', '&Tutorial'],]
        menu     = sg.Menu(menu_def, )

        #File imports
        def left(name):
            return sg.Text('{}:'.format(name), size=(10, 2))
        def mid():
            return sg.Text('Nothing imported.', size=(30, 2))
        
        repo_row       = [left('Repo Files'),   mid(), sg.FolderBrowse(key="-repo-"  )]
        combi_row      = [left('Combinations'), mid(), sg.FilesBrowse( key="-combinations-")]
        user_rules_row = [left('User Rules'),   mid(), sg.FilesBrowse( key="-user_rules-")]
        inventory_row  = [left('Inventory'),    mid(), sg.FilesBrowse( key="-inventory-") ]
        import_file    = [repo_row, combi_row, user_rules_row, inventory_row]
        
        #Main rea
        main_area = [[sg.Image(r'C:\Users\Xinling\gui\logo\final7.png'),sg.Text('HiPAD',font=('',40),text_color='#254393')],
                     [sg.Text('High Throughput Plasmid Assembly Design',font=('',12))],
                     [sg.Text('Task Name:'), sg.InputText('',key="-task_name-",size=(22,1),)],
                     [sg.Text('Assembly Method:'), sg.Drop(values=('SENAX'), key="-assembly-",size=(12,1))],
                     [sg.Frame('Import', import_file,title_color='#254393',font=(12))],
                     [sg.Text('',size=(48,2)), sg.Button('Run'), sg.Button('Cancel')]
                    ]
        
        #Overall layout
        layout = [[menu],
                  [sg.Column(main_area, element_justification='c')]
                  ]
       
        # Create the Window
        super().__init__('HiPAD', layout)
        
        #Other attributes
        self.output = None
       
    def run(self):
        window = self
        while True:
            event, self.values = window.read()
            if event == sg.WIN_CLOSED or event == 'Cancel': 
                window.close()
                break
            elif event == 'Tutorial':
                self.show_tutorial()
            elif event == 'Run':
                self.run_HiPAD()

    def __getattr__(self, attr):
        if attr == 'task_name':
            return self.values['-task_name-']
        elif attr =='assembly':
            return self.values['-assembly-']
        elif attr == 'repo':
            return self.values['-repo-']
        elif attr == 'user_rules':
            return self.values['-user_rules-']
        elif attr == 'inventory':
            return self.values['-inventory-']
        elif attr == 'combinations':
            return self.values['-combinations-']
        else:
            msg = 'HiPADGUI has no attribute "{}".'.format(attr)
            raise Exception(msg)
    
    ###############################################################################
    #Interfacing with HiPAD
    ###############################################################################
    def run_HiPAD(self):
        '''
        This method ignores inventory and user rules for now
        '''
        
        #Feedback when run button is clicked
        self.feedback()

        task_name       = self.task_name
        repo_files      = self.repo 
        combi_file      = self.combinations
        
        with open(combi_file, 'r') as file:
            data = ya.load(file, Loader=ya.FullLoader)
            combinations = data['combinations']
        
        settings = {'task_name'    : task_name,
                    'repo_files'   : repo_files,
                    'combinations' : combinations
                    }
        
        # task, settings = snx.combinatorial_synthesis(settings)
        
        settings_file = settings
        task_name, repo_files, spacers, inventory, variants, settings = snx.read_settings(settings_file)
        
        task = snx.SENAXAssembly.read_files(task_name, repo_files, spacers, inventory)

        task.assemble_variants(variants)
        task.renumber__fragments()
        
        output_folder = task.name + '_output' if type(settings_file) == dict else Path(settings_file).parent / (task.name + '_output')
        task.export_output(output_folder)
        
        self.outputs = task, settings
        
        #Retrieve folder path
        folder_path = os.getcwd() + '\\' + output_folder
        
        sg.popup_animated(None)
        
        #Feedback to inform user that download is completed and indicate folder path
        sg.popup('Download complete.',
                 'Files are saved in ' + output_folder + '.',
                 'Folder path: ' + folder_path,
                 title         = 'HiPad',
                 grab_anywhere = True,
                 )
        
    ###############################################################################
    #Other Buttons
    ###############################################################################
    def show_tutorial(self):
        tutorial = sg.popup('Repo Files: Folder containing the files with individual variant to be joined together (folder containing gb files)',
                            'Combinations: Combination of the variants to be joined together (yaml file)',
                            'User Rules: Used to specify the arrangement of variant and limit the number of instances a variant could occur in a plasmid (txt file)',
                            'Inventory: Contains all the fragments the user has in the lab. Determines whether or not PCR is required to produce the fragments required for assembly and if so, the primer sequences required.',
                            '-Note-',
                            'Repo and Combination files are compulsory.',
                            'User Rules and Inventory files are optional.',
                            title         = 'Tutorial',
                            grab_anywhere = True,
                            line_width=80,
                            )
        return tutorial
    
   
    ###############################################################################
    #Show error message and loading
    ###############################################################################
    def feedback(self):
        w = len(self.task_name)
        x = len(self.assembly)
        y = len(self.repo)
        z = len(self.combinations)

        msg = ''
        if w == 0:
            msg = msg + 'You have not given a task name.\n'
        if x == 0:
            msg = msg + 'You have not select an assembly method.\n'
        if y == 0:
            msg = msg + 'You have not import repo folder.\n'
        if z == 0:
            msg = msg + 'You have not import combinations file.\n'
        
        if w == 0 or x == 0 or y == 0 or z == 0:
            sg.popup(msg,title = 'Error')
        else:
            sg.popup_animated(sg.DEFAULT_BASE64_LOADING_GIF, no_titlebar=False, message = 'HiPad is loading, please wait for a while...',font=("Default", 12), time_between_frames=100)

        
###############################################################################
#Initialization
###############################################################################
#colour/theme
sg.LOOK_AND_FEEL_TABLE['Coolcolour'] = {'BACKGROUND': '#ACC8E2',
                                  'TEXT': 'black',
                                  'INPUT': 'white',
                                  'TEXT_INPUT': '#000000',
                                  'SCROLL': '#c7e78b',
                                  'BUTTON': ('white', '#4371A8'),
                                  'PROGRESS': ('#01826B', '#D0D0D0'),
                                  'BORDER': 1, 'SLIDER_DEPTH': 0, 'PROGRESS_DEPTH': 0,
                                    }

sg.ChangeLookAndFeel('Coolcolour')

if __name__ == '__main__':
    window = HiPADGUI()
    window.run()