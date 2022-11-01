const { contextBridge, ipcRenderer } = require('electron');

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