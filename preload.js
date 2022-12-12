const { contextBridge, ipcRenderer } = require('electron');
const fs = require("fs");
const path = require('path');

window.addEventListener('DOMContentLoaded', () => {
    const replaceText = (selector, text) => {
        const element = document.getElementById(selector)
        if (element) element.innerText = text
    }
  
    for (const type of ['crns']) {
        replaceText(`${type}-version`, fs.readFileSync(path.join(__dirname, 'VERSION'), 'utf8'));
    }
    const doing = document.getElementById('doing');
    const progress = document.getElementById('progress');
    ipcRenderer.on('update_progress', (_event, value) => {
        progress.innerText = value;
    });
    ipcRenderer.on('update_doing', (_event, value) => {
        doing.innerText = value;
    });
});

contextBridge.exposeInMainWorld(
    'electron',
    {
        open: () => ipcRenderer.invoke('OPEN_FILE').then((result) => {
            if(result.canceled == false) {
                var paths = result.filePaths;
                console.log(paths);
                return paths[0];
            } else {
                return undefined;
            }
            
        }),
        save: () => ipcRenderer.invoke('SAVE_FILE').then((result) => {
            if(result.canceled == false) {
                var path = result.filePath;
                console.log(result)
                return path;
            } else {
                return undefined;
            }
        })

    }
)