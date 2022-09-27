# -*- mode: python ; coding: utf-8 -*-
import site

_p = site.getsitepackages()[0]
if "site-packages" in _p:
    _p = os.path.join(site.getsitepackages()[0],'eel','eel.js')
else:
    _P = os.path.join(site.getsitepackages()[0], 'lib', 'site-packages','eel','eel.js')

block_cipher = None


a = Analysis(
    ['main.py'],
    pathex=[],
    binaries=[],
    datas=[(_p, 'eel'), ('static', 'static')],
    hiddenimports=['bottle_websocket'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [('u', None, 'OPTION')],
    exclude_binaries=True,
    name='main',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='main',
)
