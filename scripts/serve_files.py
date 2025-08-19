#!/usr/bin/env python3
"""
Simple HTTP server to serve AB1 files for the browser application
Run this from the project root directory: python3 scripts/serve_files.py
"""

import http.server
import socketserver
import os
import sys
from pathlib import Path

# Set the port
PORT = 8000

# Get the project root directory
PROJECT_ROOT = Path(__file__).parent.parent
os.chdir(PROJECT_ROOT)

class CustomHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # Add CORS headers
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', '*')
        super().end_headers()

    def do_OPTIONS(self):
        self.send_response(200)
        self.end_headers()

    def guess_type(self, path):
        """Override to handle .ab1 files properly"""
        if str(path).endswith('.ab1'):
            return 'application/octet-stream'
        return super().guess_type(path)

if __name__ == "__main__":
    print(f"🚀 Starting HTTP server on port {PORT}")
    print(f"📂 Serving files from: {PROJECT_ROOT}")
    print(f"🌐 Browser app: http://localhost:{PORT}/browser/")
    print(f"📁 AB1 files: http://localhost:{PORT}/input/")
    print("Press Ctrl+C to stop")
    
    with socketserver.TCPServer(("", PORT), CustomHTTPRequestHandler) as httpd:
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\n🛑 Server stopped")
            sys.exit(0)
