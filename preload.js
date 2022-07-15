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
            
        })
    }
)