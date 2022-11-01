const { app, BrowserWindow, dialog, ipcMain } = require('electron');
const proc = require('child_process');
const url = require('url');
const path = require('path');
const isDev = require('electron-is-dev');
//var icmp = require('icmp');

command = path.join(__dirname, 'dist', 'main', 'main')
args = []
if(isDev) {
    console.log("DEV-Build");
    command = 'python';
    args = ['-u','main.py']
} 

//let child = proc.spawn('python', ['-u','main.py']);
let child = proc.spawn(command, args);

let eel_started = false;

child.stdout.on('data', function(data) {
    console.log('stdout: ' + data);

    if (data.toString().includes('[eel]: Start')) {
        eel_started = true;
    }
});

child.stderr.on('data', function(data) {
    console.log('stderr: ' + data);
});

async function createWindow () {
    const win = new BrowserWindow({
        width: 800,
        height: 600,
        webPreferences: {
            enableRemoteModule: true,
            preload: path.join(app.getAppPath(), 'preload.js')
        }
    });
    while(eel_started == false) {
        await sleep(500);
        win.loadURL('http://localhost:8000/');
    }
    win.loadURL('http://localhost:8000/'); //.catch();
}

// icmp.send('http://localhost:8000/', "Hey, I'm sending a message!")
//     .then(obj => {
//         console.log(obj.open ? 'Done' : 'Failed')
//     })
//     .catch(err => console.log(err));

function sleep(ms) {
    return new Promise((resolve) => {
      setTimeout(resolve, ms);
    });
}

app.whenReady().then(() => {
    createWindow()

    app.on('activate', () => {
        if (BrowserWindow.getAllWindows().length === 0) {
            createWindow()
        }
    })
})

app.on('window-all-closed', () => {
    child.kill('SIGINT')
    app.quit()
})


ipcMain.handle('OPEN_FILE', async (event) => {
    const filePath = await dialog.showOpenDialog({filters: [
        { name: 'Network', extensions: ['txt', 'xml', 'sbml'] }
    ]});
    return filePath;
});

ipcMain.handle('SAVE_FILE', async (event) => {
    const filePath = await dialog.showSaveDialog();
    return filePath;
});